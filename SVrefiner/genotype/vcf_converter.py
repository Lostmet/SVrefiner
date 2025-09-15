import os
import ast
import pysam
from typing import NamedTuple
import pandas as pd
from ..utils.logging_utils import get_logger


logger = get_logger(__name__)


class PlinkFiles(NamedTuple):
    bed: str
    bim: str
    fam: str

def nrSV_vcf_generate(nrSV_csv: str, config, nrSV_list, nSV_name: str):
    '''
    Generate a new VCF-style table from nrSV metadata and original VCF genotype info.
    '''
    df = pd.read_csv(nrSV_csv) 
    nrSV_count = df.shape[0]
    output_cache = {}  # Dictionary to store SV genotype results

    # Get the output VCF file name
    output_vcf = os.path.join(config.output_dir, nSV_name)

    # Read the original VCF file
    with pysam.VariantFile(config.vcf_file, "r") as vcf_in:
        for item in nrSV_list:
            parts = item['group'].split('_')
            chrom = parts[1]
            pos = int(parts[3])
            sv_indexes = item['SV']

            # Initialize genotype storage for this group
            output_cache[item['group']] = {}

            # Traverse all SVs at this VCF position
            sv_records = []
            for record in vcf_in.fetch(chrom, pos - 1, pos):  # VCF is 1-based, fetch is 0-based
                if record.pos == pos:
                    sv_records.append(record)

            # Extract genotype information for given SV indexes
            for sv_idx in sv_indexes:
                if sv_idx <= len(sv_records):  # Ensure the index is valid
                    sv_record = sv_records[sv_idx - 1]  # SV index starts from 1; Python list is 0-based

                    # Get GT (genotype) info for all samples
                    gt_data = {sample: sv_record.samples[sample]['GT'] for sample in sv_record.samples}

                    # Store GT info into cache
                    output_cache[item['group']][f"SV{sv_idx}"] = gt_data 

    # Retrieve sample names from the first available entry in output_cache
    sample_names = []
    if output_cache:
        first_group = next(iter(output_cache))  # First group name
        first_sv = next(iter(output_cache[first_group])) if output_cache[first_group] else None  # First SV_x

        if first_sv:
            sample_names = list(output_cache[first_group][first_sv].keys())  # Extract sample name list

    # Write output VCF-style rows
    with open(output_vcf, "a") as f_out:
        for _, row in df.iterrows():
            group_name_parts = row['chromosome_group'].split('_')
            group_name = '_'.join(group_name_parts[:-1])
            nrSV_id = row['chromosome_group']
            chrom = group_name_parts[1]
            meta = ast.literal_eval(row['meta_array'])
            sv_id = int(row['sequence_id'][3:])
            ref = meta['ref']
            alt = meta['alt']
            pos = meta['pos']
            try:
                gt_info = output_cache[group_name][f'SV{sv_id}']
            except KeyError as e:
                logger.warning(f"KeyError: {e} not found in {group_name} when processing nrSV")
                gt_info = {}
            gt_values = ["/".join(map(str, gt_info.get(sample, ('.', '.')))) for sample in sample_names]
            f_out.write(f"{chrom}\t{pos}\t{nrSV_id}\t{ref}\t{alt}\t.\t.\ttype=nrSV\tGT\t" + "\t".join(gt_values) + "\n")

    return nrSV_count
