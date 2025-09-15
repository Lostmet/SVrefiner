# === Standard library ===
import os
import re
import sys
import csv
import time

# === Third-party libraries ===
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

# === Project modules ===
from ..utils.logging_utils import get_logger


logger = get_logger(__name__)


def sample_name_contract(vcf_name: str):
    """Extract sample names from VCF file header."""
    vcf_file = pysam.VariantFile(vcf_name, 'r')
    sample_names = list(vcf_file.header.samples)
    return sample_names


def _dfile_sort_key_true_style(fname: str):
    """
    Keep consistent with load_encoded_gt_from_matrix_dir:
    only use the 2nd and 3rd segments as numeric sort keys.
    Example: Group_1_23_456_D_matrix.csv -> key = (1, 23)
    """
    parts = fname.split('_')
    try:
        return (int(parts[1]), int(parts[2]))
    except Exception:
        return (1 << 60, 1 << 60)


def process_vcf_to_x_matrix(vcf_dir: str, output_dir: str, writedown: bool = False):
    """
    Extract GT matrix from VCF.
    writedown=False: do not write X, compute T in memory and encode GT, return gt_buffer.
    writedown=True : keep old logic (write X to disk, gt_buffer returns None).
    """

    vcf_file = pysam.VariantFile(os.path.join(vcf_dir, "oSV.vcf"), 'r')
    sample_names = list(vcf_file.header.samples)

    # === Scan VCF records (variants) ===
    gt_data = []
    with tqdm(desc="VCF→GT (scan variants)", unit="var", mininterval=0.1) as pbar:
        for record in vcf_file:
            gt_row = [record.contig, record.pos]
            for s in sample_names:
                gt = record.samples[s].get("GT")
                if gt is None or len(gt) != 2 or any(a is None for a in gt):
                    gt_row.append(-999)
                else:
                    a, b = gt
                    if   (a, b) == (0, 0): gt_row.append(0)
                    elif (a, b) in ((0, 1), (1, 0)): gt_row.append(1)
                    elif (a, b) == (1, 1): gt_row.append(2)
                    else: gt_row.append(-1)
            gt_data.append(gt_row)
            pbar.update(1)

    header = ["#CHROM", "POS"] + sample_names
    df_vcf = pd.DataFrame(gt_data, columns=header)

    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    csv_files = [f for f in os.listdir(output_dir) if re.match(r'Group_\d+_\d+_\d+_D_matrix\.csv', f)]
    # Consistent with True style sorting (ensure _dfile_sort_key_true_style is defined)
    csv_files = sorted(csv_files, key=_dfile_sort_key_true_style)

    df_vcf["#CHROM"] = df_vcf["#CHROM"].astype(str)

    gt_buffer = [] if not writedown else None

    # === Process groups (aligned to D, one group at a time) ===
    total_groups = len(csv_files)
    with tqdm(total=total_groups, desc="Genotyping rSV", unit="group", mininterval=0.2) as pbar_groups:
        for csv_file in csv_files:
            m = re.match(r'Group_(\d+)_(\d+)_(\d+)_D_matrix\.csv', csv_file)
            if not m:
                pbar_groups.update(1)
                continue
            chrom, number, pos = m.groups()
            chrom, pos = str(chrom), int(pos)

            vcf_row_index = df_vcf[(df_vcf["#CHROM"] == chrom) & (df_vcf["POS"] == pos)].index
            if vcf_row_index.empty:
                logger.warning(f"Warning: #CHROM {chrom}, POS {pos} not found in VCF for {csv_file}")
                pbar_groups.update(1)
                continue

            start_idx = vcf_row_index[0]
            d_path = os.path.join(output_dir, csv_file)
            df_d_full = pd.read_csv(d_path)  # first column: block key; subsequent columns: rSV IDs
            num_rows_to_extract = len(df_d_full)

            vcf_subset = df_vcf.iloc[start_idx: start_idx + num_rows_to_extract].reset_index(drop=True)
            gt_matrix = vcf_subset[sample_names]
            gt_np = gt_matrix.to_numpy(dtype=np.int64, copy=False)

            if writedown:
                # Old logic: write X to disk (unchanged)
                output_file = os.path.join(output_dir, csv_file.replace("_D_matrix.csv", "_X_matrix.csv"))
                header_row = [f"{chrom}_{pos}"] + sample_names
                first_col = df_d_full.iloc[:, 0].astype(str).tolist()
                with open(output_file, "w", newline="", buffering=1024*1024) as f:
                    w = csv.writer(f, lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
                    w.writerow(header_row)
                    for r in range(gt_np.shape[0]):
                        row = [first_col[r]]
                        row.extend(gt_np[r, :].tolist())
                        w.writerow(row)
            else:
                # New logic: do not write to disk, directly compute T and encode → append to gt_buffer
                d_data = df_d_full.iloc[:, 1:].to_numpy(dtype=np.int64, copy=False)  # (rows, m)
                t_mat  = d_data.T.dot(gt_np)                                        # (m, n_samples)

                enc = np.full(t_mat.shape, "./.", dtype=object)
                enc[t_mat == 0] = "0/0"
                enc[t_mat == 1] = "1/0"
                enc[t_mat >= 2] = "1/1"

                gt_buffer.extend(enc.tolist())

            pbar_groups.update(1)

    return sample_names, gt_buffer
