import os
import pandas as pd
import re
import pysam
import numpy as np
from ..utils.logging_utils import get_logger


logger = get_logger(__name__)


def sample_name_contract(vcf_name: str):
    """Extract sample names from VCF file header."""
    vcf_file = pysam.VariantFile(vcf_name, 'r')
    sample_names = list(vcf_file.header.samples)
    return sample_names


def process_vcf_to_x_matrix(vcf_dir: str, output_dir: str):
    """Extract GT matrix from VCF and generate X_matrix files by matching D_matrix."""
    vcf_file = pysam.VariantFile(os.path.join(vcf_dir, "oSV.vcf"), 'r')
    sample_names = list(vcf_file.header.samples)

    gt_data = []
    for record in vcf_file:
        gt_row = [record.contig, record.pos]
        for sample in sample_names:
            gt = record.samples[sample]["GT"]
            if gt is None or len(gt) != 2 or any(allele is None for allele in gt):
                gt_row.append("./.")
            else:
                gt_row.append(f"{gt[0]}/{gt[1]}")
        gt_data.append(gt_row)

    header = ["#CHROM", "POS"] + sample_names
    df_vcf = pd.DataFrame(gt_data, columns=header)

    def transform_gt(gt):
        if gt == './.': return -999
        elif gt == '0/0': return 0
        elif gt == '1/0' or gt == '0/1': return 1
        elif gt == '1/1': return 2
        else: return -1

    for sample in sample_names:
        df_vcf[sample] = df_vcf[sample].apply(transform_gt)

    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    csv_files = [f for f in os.listdir(output_dir) if re.match(r'Group_\d+_\d+_\d+_D_matrix\.csv', f)]

    for csv_file in csv_files:
        match = re.match(r'Group_(\d+)_(\d+)_(\d+)_D_matrix\.csv', csv_file)
        if not match:
            continue

        chrom, number, pos = match.groups()
        pos = int(pos)
        chrom = str(chrom)

        vcf_row_index = df_vcf[(df_vcf["#CHROM"].astype(str) == chrom) & (df_vcf["POS"] == pos)].index

        if vcf_row_index.empty:
            logger.warning(f"Warning: #CHROM {chrom}, POS {pos} not found in VCF for {csv_file}")
            continue

        start_idx = vcf_row_index[0]
        csv_data = pd.read_csv(os.path.join(output_dir, csv_file), usecols=[0])
        num_rows_to_extract = len(csv_data)

        vcf_subset = df_vcf.iloc[start_idx: start_idx + num_rows_to_extract].reset_index(drop=True)
        gt_matrix = vcf_subset[sample_names]

        csv_data.columns = [f"{chrom}_{pos}"]
        csv_data = pd.concat([csv_data, gt_matrix], axis=1)

        output_file = os.path.join(output_dir, csv_file.replace("_D_matrix.csv", "_X_matrix.csv"))
        csv_data.to_csv(output_file, index=False)

        if not os.path.exists(output_file):
            logger.warning(f"Warning: {output_file} not created!")
        else:
            continue

    return sample_names


def compute_t_matrix(output_dir: str):
    """Compute T_matrix = D × X and save to disk."""
    output_dir = os.path.abspath(output_dir)
    files = os.listdir(output_dir)

    d_files = [f for f in files if re.match(r'Group_\d+_\d+_\d+_D_matrix\.csv', f)]
    x_files = [f for f in files if re.match(r'Group_\d+_\d+_\d+_X_matrix\.csv', f)]

    d_files.sort()
    x_files.sort()

    if len(d_files) != len(x_files):
        logger.warning("Warning: D_matrix and X_matrix file counts do not match.")
        return

    for d_file, x_file in zip(d_files, x_files):
        match_d = re.match(r'Group_(\d+)_(\d+)_(\d+)_D_matrix\.csv', d_file)
        match_x = re.match(r'Group_(\d+)_(\d+)_(\d+)_X_matrix\.csv', x_file)

        if not match_d or not match_x:
            logger.warning(f"Warning: Skipping unmatched files: {d_file}, {x_file}")
            continue

        df_d = pd.read_csv(os.path.join(output_dir, d_file))
        df_x = pd.read_csv(os.path.join(output_dir, x_file))

        d_data = df_d.iloc[:, 1:].values
        x_data = df_x.iloc[:, 1:].values

        t_matrix = np.dot(d_data.T, x_data)
        t_df = pd.DataFrame(t_matrix, index=df_d.columns[1:], columns=df_x.columns[1:])
        first_column_name = df_x.columns[0]
        t_df.index.name = first_column_name

        t_matrix_file = os.path.join(output_dir, d_file.replace("_D_matrix.csv", "_T_matrix.csv"))
        t_df.to_csv(t_matrix_file, index=True, header=True)
