import pysam
from tqdm import tqdm
import concurrent.futures
from ..config import Config
from typing import Dict, List, Tuple
from ..utils.logging_utils import get_logger
from .vcf_parser import VariantGroup, Variant

logger = get_logger(__name__)


def get_max_insertions(variants: List[Variant], start: int) -> Tuple[Dict[int, int], bool]:
    """Identify maximum insertion length per position and detect overlapping insertions."""
    max_insertions = {}
    pos_set = []
    has_insertion = False

    for variant in variants:
        rel_start = variant.start - start
        if len(variant.ref) == 1 and len(variant.alt[0]) > 1:
            max_insertions[rel_start] = max(max_insertions.get(rel_start, 0), len(variant.alt[0]))
            pos_set.append(variant.start)

    if len(pos_set) != len(set(pos_set)):
        has_insertion = True

    return max_insertions, has_insertion


def adjust_reference_for_insertions(ref_seq: str, max_insertions: Dict[int, int]) -> str:
    """Insert dashes into reference sequence to align with maximum insertions."""
    ref_list = list(ref_seq)
    ins = 0
    for pos, max_ins_length in sorted(max_insertions.items()):
        blanks = list("-" * (max_ins_length - 1))
        for i, char in enumerate(blanks):
            ref_list.insert(pos + ins + i, char)
        ins += max_ins_length - 1
    return "".join(ref_list)


def adjust_variants_for_insertions(
    ref_seq: str, max_insertions: Dict[int, int], variant: Variant, start: int
) -> Tuple[str, Dict[str, int]]:
    ref_list = list(ref_seq)
    ins = 0
    rel_start = variant.start - start
    rel_end = variant.end - start + 1
    poly_ins = {}

    if len(variant.ref) > 1 and len(variant.alt[0]) == 1:
        ref_list[rel_start:rel_end - 1] = list((len(variant.ref) - 1) * "-")
        for pos, max_ins_length in sorted(max_insertions.items()):
            blanks = list("-" * (max_ins_length - 1))
            for i, char in enumerate(blanks):
                ref_list.insert(pos + ins + i, char)
            ins += max_ins_length - 1

    elif len(variant.ref) == 1 and len(variant.alt[0]) > 1:
        for pos, max_ins_length in sorted(max_insertions.items()):
            if rel_start != pos:
                blanks = list("-" * (max_ins_length - 1))
                for i, char in enumerate(blanks):
                    ref_list.insert(pos + ins + i, char)
                ins += max_ins_length - 1
            else:
                alt = list(variant.alt[0][1:])
                blanks = list("-" * (max_ins_length - len(variant.alt[0])))
                position = 0
                for char in alt:
                    ref_list.insert(pos + ins + position, char)
                    position += 1
                for char in blanks:
                    ref_list.insert(pos + ins + position, char)
                    position += 1
                if blanks:
                    poly_ins["start"] = pos + ins
                    poly_ins["end"] = pos + ins + position
                    poly_ins["pos"] = start + 1
                ins += max_ins_length - 1

    else:
        ref_list[rel_start:rel_end - 1] = list(variant.alt)
        for pos, max_ins_length in sorted(max_insertions.items()):
            blanks = list("-" * (max_ins_length - 1))
            for i, char in enumerate(blanks):
                ref_list.insert(pos + ins + i, char)
            ins += max_ins_length - 1

    return "".join(ref_list), poly_ins


def generate_fasta_sequences(
    config: Config, variant_groups: Dict[str, List[VariantGroup]], total_groups: int
) -> Tuple[str, Dict[str, bool], List[Dict[str, int]]]:
    """Generate pre-aligned FASTA sequences and collect insertion information."""
    output_fasta = f"{config.output_dir}/variants_pre_aligned.fasta"
    ref_genome = pysam.FastaFile(config.ref_fasta)
    has_insertion_dict = {}
    poly_ins_list = []

    def process_group(chrom, i, group):
        start = group.start - 1
        end = group.end
        ref_seq = ref_genome.fetch(chrom, start, end).upper()
        max_insertions, has_insertion = get_max_insertions(group.variants, start)
        has_insertion_dict[f"Group_{chrom}_{i}_{group.start}"] = has_insertion
        ref_seq_adjusted = adjust_reference_for_insertions(ref_seq, max_insertions)
        seq_length = len(ref_seq_adjusted)
        variant_results = []

        for variant in group.variants:
            var_seq, poly_ins = adjust_variants_for_insertions(ref_seq, max_insertions, variant, start)
            var_id = f"Variant_{chrom}_{i}_{variant.start}_{variant.end}"
            variant_results.append((var_id, var_seq, poly_ins))

        bp_to_process = (len(group.variants) + 1) * seq_length
        return chrom, i, group, ref_seq_adjusted, variant_results, bp_to_process

    try:
        with open(output_fasta, 'w') as fasta_out:
            with tqdm(total=total_groups, desc="Pre-aligning variant groups", unit="group") as pbar:
                futures = []
                with concurrent.futures.ThreadPoolExecutor(max_workers=config.threads) as executor:
                    for chrom, groups in variant_groups.items():
                        for i, group in enumerate(groups, 1):
                            futures.append(executor.submit(process_group, chrom, i, group))

                    seen = set()
                    group_bp_all = 0
                    for future in concurrent.futures.as_completed(futures):
                        chrom, i, group, ref_seq_adjusted, variant_results, group_bp = future.result()
                        group_bp_all += group_bp

                        fasta_out.write(f">Group_{chrom}_{i}_{group.start}\n{ref_seq_adjusted}\n")

                        for var_id, var_seq, poly_ins in variant_results:
                            fasta_out.write(f">{var_id}\n{var_seq}\n")
                            if poly_ins:
                                dict_frozenset = frozenset(poly_ins.items())
                                if dict_frozenset not in seen:
                                    poly_ins_list.append(poly_ins)
                                    seen.add(dict_frozenset)

                        pbar.update(1)

        return output_fasta, has_insertion_dict, poly_ins_list

    finally:
        ref_genome.close()
