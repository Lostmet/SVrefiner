import os
import pysam
import click
from tqdm import tqdm
from dataclasses import dataclass
from ..exceptions import InputError
from typing import Dict, List, Tuple
from ..utils.logging_utils import get_logger, log_tqdm_summary


logger = get_logger(__name__)


@dataclass
class Variant:
    chrom: str
    start: int
    end: int
    ref: str
    alt: List[str]

class VariantGroup:
    def __init__(self, chrom: str):
        self.chrom = chrom
        self.variants = []
        self.start = float('inf')
        self.end = 0
        
    def add_variant(self, variant: Variant):
        self.variants.append(variant)
        self.start = min(self.start, variant.start)
        self.end = max(self.end, variant.end)
        
    def overlaps(self, variant: Variant) -> bool:
        return (variant.start <= self.end and 
                variant.end >= self.start)


def parse_chrom(chrom: str) -> int:
    """
    Convert the chromosome string to an integer that reflects a more natural
    ordering (e.g., 1 < 2 < ... < 10 < X < Y, etc.). 
    """
    # Remove common prefix like 'chr' if present
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    
    try:
        return int(chrom)
    except ValueError:
        logger.error(f"Non-standard chromosome name detected, chrom: {chrom}. Please encode it as an integer manually.")
        raise click.Abort()

def process_variants(config) -> List[VariantGroup]:
    """
    Process a VCF file, sorting all variants by chromosome order 
    (using parse_chrom) and by their start coordinate, and then grouping 
    overlapping variants into VariantGroup objects.
    """
    try:
        vcf = pysam.VariantFile(config.vcf_file)
    except Exception as e:
        raise InputError(f"Failed to open VCF file: {str(e)}")
        
    all_variants = []
    single_group = [] 
    inv_group = [] 
    var_not_align_count = 0 
    var_not_align_list = []
    var_not_standard_count = 0
    snp_group = []
    var_not_aligned_group = []

    base = {"A","T","C","G","N"}
    def is_simple_base(seq: str) -> bool:
        return set(seq.upper()) <= base

    with tqdm(desc="Reading variants",unit='variants') as pbar:
        try:
            var_bp_all = 0
            var_bp_max = 0
            for record in vcf.fetch():
                if list(record.alts) == ['<INV>']:
                    inv = Variant(
                        chrom=record.chrom,
                        start=record.pos,
                        end=record.pos + len(record.ref) - 1,
                        ref=record.ref,
                        alt=list(record.alts)
                    )
                    inv_group.append(inv) 
                else:
                    if_snp = False
                    if_not_aligned = False
                    # Modify pos，ref and alt, left-aligned
                    chrom = record.chrom
                    start = record.pos
                    pos = start
                    ref = record.ref
                    alt = list(record.alts)
                    ref_len = len(ref)
                    alt_len = len(alt[0])
                    end = pos + ref_len -1
                    if not is_simple_base(ref) or not is_simple_base(alt[0]):
                        if var_not_standard_count == 0:
                            var_not_align_list.append(f'chrom:{chrom}, pos:{pos}')
                            logger.warning(
                                f"Detected illegal REF/ALT base at chrom:{chrom}, pos:{pos}. "
                                f"Possible wrong VCF. Please check."
                            )
                        else:
                            var_not_align_list.append(f'chrom:{chrom}, pos:{pos}')
                        var_not_standard_count += 1
                        var_not_aligned = Variant(
                            chrom=chrom,
                            start=pos,
                            end=end,
                            ref=ref,
                            alt=alt
                        )
                        var_not_aligned_group.append(var_not_aligned)
                        if_not_aligned = True

                    elif ref_len == 1 and alt_len == 1: # SNP
                        snp = Variant(
                            chrom=chrom,
                            start=pos,
                            end=end,
                            ref=ref,
                            alt=alt
                        )
                        snp_group.append(snp) # to new group 
                        if_snp = True                     

                    elif ref_len != 1 and alt_len != 1: # not left-aligned or complex variant
                        if var_not_align_count == 0:
                            var_not_align_list.append(f'chrom:{chrom}, pos:{pos}')
                            logger.warning(f"Warning: Unaligned variant detected, chrom: {chrom}, pos: {pos}. \
                                           \nPlease check; the software will save these variants into nSV.vcf")
                        else:
                            var_not_align_list.append(f'chrom:{chrom}, pos:{pos}')

                        var_not_align_count += 1
                        var_not_aligned = Variant(
                            chrom=chrom,
                            start=pos,
                            end=end,
                            ref=ref,
                            alt=alt
                        )
                        var_not_aligned_group.append(var_not_aligned)
                        if_not_aligned = True

                    elif ref[0] != alt[0][0]:
                        if var_not_standard_count == 0:
                            var_not_align_list.append(f'chrom:{chrom}, pos:{pos}')
                            logger.warning(f"Detected variant where either REF or ALT is 1 bp long, but the first base differs, indicating possible complex substitution events\
                                            \nchrom: {chrom}, pos: {pos}. Please check the VCF file.")
                        else:
                            var_not_align_list.append(f'chrom:{chrom}, pos:{pos}')
                        var_not_standard_count += 1
                        var_not_aligned = Variant(
                            chrom=chrom,
                            start=pos,
                            end=end,
                            ref=ref,
                            alt=alt
                        )
                        var_not_aligned_group.append(var_not_aligned)
                        if_not_aligned = True
                    
                    
                    if (if_snp == False) and (if_not_aligned == False):
                        variant = Variant(
                            chrom=chrom,
                            start=pos,
                            end=end,
                            ref=ref,
                            alt=alt
                        )
                        var_bp = abs(len(variant.ref)-len(variant.alt[0]))
                        var_bp_max = max(var_bp_max, var_bp)
                        var_bp_all += var_bp
                        all_variants.append(variant)
                pbar.update(1)
            pbar.close()
            log_tqdm_summary(pbar, logger)
        except Exception as e:
            raise InputError(f"Error reading VCF file: {str(e)}")
        finally:
            vcf.close()

    all_variants.sort(key=lambda v: (parse_chrom(v.chrom), v.start))

    # Filter abnormal variants
    variant_groups = []
    inv_groups = []
    snp_groups = []
    var_not_aligned_groups = []
    current_group = None
    inv_count = len(inv_group)
    snp_count = len(snp_group)
    var_abnormal_count = var_not_align_count +  var_not_standard_count
    variant_count = len(all_variants) + len(var_not_aligned_group) + inv_count + snp_count

    if var_abnormal_count != 0:
        log_path = os.path.join(config.output_dir, "VCF_normalization_failed.log", )
        logger.warning(
    f"Warning: A total of {var_not_align_count + var_not_standard_count} variants may not be properly left-aligned.\n"
    "Please refer to the README for normalization instructions and double-check your VCF input.\n"
    f"(Detailed log saved at: {log_path})"
)
        with open(log_path, "w") as f:
            f.write(f"{var_not_align_count} variants have redundant REF and ALT sequence overlap (i.e., both REF and ALT are longer than 1 bp)\n")
            f.write(f"{var_not_standard_count} variants have non-matching bases at the first position of REF and ALT\n")
            f.write("These issues may indicate that your VCF is not left-aligned or not normalized properly.\n")
            f.write(f"Abnormal variants account for {100*(var_not_align_count + var_not_standard_count)/variant_count:.2f}% of all detected variants.\n")
            f.write("Please check the following variant records:\n")
            for item in var_not_align_list:
                f.write(item + "\n")

    # Process inv_group
    if inv_group:
        for inv in inv_group:
            inv_current_group = VariantGroup(inv.chrom)
            inv_current_group.add_variant(inv)
            inv_groups.append(inv_current_group)

        logger.info(f"Excluding {inv_count} inversion(s)")

    # Process snp_group
    if snp_group:
        for snp in snp_group:
            snp_current_group = VariantGroup(snp.chrom)
            snp_current_group.add_variant(snp)
            snp_groups.append(snp_current_group)

        logger.warning(f"Excluding {snp_count} SNP(s)") 

    # Process var_not_aligned_group
    if var_not_aligned_group:
        for var in var_not_aligned_group:
            var_current_group = VariantGroup(var.chrom)
            var_current_group.add_variant(var)
            var_not_aligned_groups.append(var_current_group)

        logger.warning(f"Excluding {var_abnormal_count} variant(s) of not normalized")

    # Grouping variants
    for variant in all_variants:
        if current_group is None:
            # First variant
            current_group = VariantGroup(variant.chrom)
            current_group.add_variant(variant)
        else:
            # If chromosome is different or they do not overlap, start a new group
            if (current_group.chrom != variant.chrom or 
                not current_group.overlaps(variant)):
                variant_groups.append(current_group)
                current_group = VariantGroup(variant.chrom)
                current_group.add_variant(variant)
            else:
                # Still in the same group
                current_group.add_variant(variant)

        # dummy update (replace tqdm update)
        pass

    # Add the last group if it exists
    if current_group is not None:
        variant_groups.append(current_group)


    multi_group = []
    multi_var_bp = 0
    multi_var_bp_max = 0
    single_sv_count = 0
    multi_var_bp_all = 0
    variant_max = None

    # inv_groups to single
    for i in range(len(inv_groups)):
        single_group.append(inv_groups[i])
    single_sv_count += len(inv_groups)

    for i in range(len(var_not_aligned_groups)):
        single_group.append(var_not_aligned_groups[i])
    single_sv_count += len(var_not_aligned_groups)
    
    ### multi_group1:{'chrom': '2', 'variants': [Variant(chrom='2', start=906670, end=906670, ref='A', alt=['ATATATATATATA'], samples={'SL001_SL001':
    for i in range(len(variant_groups)):
        if len(variant_groups[i].variants) == 1:
            single_group.append(variant_groups[i])
            single_sv_count += 1
        else:
            multi_group.append(variant_groups[i])
            for n in range(len(variant_groups[i].variants)):
                variant = variant_groups[i].variants[n]
                ref_bp = len(variant.ref)
                alt_bp = len(variant.alt[0])
                multi_var_bp = abs(ref_bp - alt_bp)
                multi_var_bp_all += multi_var_bp
                if multi_var_bp_max < multi_var_bp:
                    multi_var_bp_max = multi_var_bp
                    variant_max = variant

    percentage_sv_overlapped = (1- single_sv_count/variant_count)*100 # get overlapped SV percentage (exclude SNP)
    single_group.sort(key=lambda v: (parse_chrom(v.chrom), v.start))
    logger.info(f"Total base pairs to be processed: {var_bp_all:,}, max per variant: {var_bp_max:,}")
    logger.info(f"After grouping: {multi_var_bp_all:,} base pairs, max per variant: {multi_var_bp_max:,}")
    if variant_max:
        logger.info(f"Max variant located at: chromosome {variant_max.chrom}, position {variant_max.start:,}")
    logger.info(f"{single_sv_count:,} SV(s) excluded due to lack of overlap")
    logger.info(f"Percentage of overlapping SVs: {percentage_sv_overlapped:.2f}%")

    if variant_max == None:
        click.echo("Have no overlapped SV")
        raise click.Abort()

    return multi_group, single_sv_count, multi_var_bp_all, percentage_sv_overlapped, single_group, inv_count, variant_count, snp_groups


def filter_vcf(config, single_group, snp_groups):
    nvcf_name = "nSV.vcf"
    if single_group:
        output_vcf = os.path.join(config.output_dir, nvcf_name)
        with pysam.VariantFile(config.vcf_file, "r") as vcf_in:
            header_text = str(vcf_in.header)
        
        with open(output_vcf, "w") as f_out:
            f_out.write(header_text)
            
            with pysam.VariantFile(config.vcf_file, "r") as vcf_in:
                for i in tqdm(range(len(single_group)), desc="Generating nSV: ", unit='variants'):
                    variant = single_group[i].variants[0]
                    chrom, pos, ref, alt = variant.chrom, variant.start, variant.ref, variant.alt[0]
                    
                    for rec in vcf_in.fetch(chrom, pos - 1, pos + 1):
                        if rec.pos == pos and rec.ref == ref and rec.alts and rec.alts[0] == alt:
                            vcf_line = str(rec).strip()
                            f_out.write(vcf_line + '\n')
                            break

        logger.info(f"Non-overlapping SVs(nSVs) saved at {output_vcf}")
    if snp_groups:
        snpvcf_name = "SNP.vcf"
        output_vcf = os.path.join(config.output_dir, snpvcf_name)
        with pysam.VariantFile(config.vcf_file, "r") as vcf_in:
            header_text = str(vcf_in.header)
        
        with open(output_vcf, "w") as f_out:
            f_out.write(header_text)
            
            with pysam.VariantFile(config.vcf_file, "r") as vcf_in:
                for i in tqdm(range(len(snp_groups)), desc="Generating SNP: ", unit='variants'):
                    variant = snp_groups[i].variants[0]
                    chrom, pos = variant.chrom, variant.start
                    
                    for rec in vcf_in.fetch(chrom, pos - 1, pos + 1):
                        if rec.pos == pos:
                            vcf_line = str(rec).strip()
                            f_out.write(vcf_line + '\n')
                            break        
        logger.info(f"SNP VCF file saved at {output_vcf}")

    return nvcf_name


def make_oSV(config, multi_group):
    ovcf_name = "oSV.vcf"
    if multi_group:
        output_vcf = os.path.join(config.output_dir, ovcf_name)
        with pysam.VariantFile(config.vcf_file, "r") as vcf_in:
            header_text = str(vcf_in.header)
        
        with open(output_vcf, "w") as f_out:
            f_out.write(header_text)
            
            with pysam.VariantFile(config.vcf_file, "r") as vcf_in:
                for i in tqdm(range(len(multi_group)), desc="Generating oSV: ", unit='groups'):
                    variant_group = multi_group[i].variants
                    for variant in variant_group:
                        chrom, pos, ref, alt = variant.chrom, variant.start, variant.ref, variant.alt[0]
                        
                        for rec in vcf_in.fetch(chrom, pos - 1, pos + 1):
                            if rec.pos == pos and rec.ref == ref and rec.alts and rec.alts[0] == alt:
                                vcf_line = str(rec).strip()
                                f_out.write(vcf_line + '\n')
                                break

        logger.info(f"Overlapping SVs(oSVs) saved at {output_vcf}")    
