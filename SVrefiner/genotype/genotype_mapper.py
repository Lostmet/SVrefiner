# === Standard library ===
import os
import re
import sys
import ast
import csv
import time
import logging
from pathlib import Path
from typing import List, Tuple, Optional

# === Third-party libraries ===
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

# === Project modules ===
from ..utils.logging_utils import get_logger


logger = get_logger(__name__)
csv.field_size_limit(sys.maxsize)


def yes_no_to_bool(value):
    return value.strip().lower() in ['yes', 'y', 'true', '1']


def compute_t_matrix(output_dir: str):
    """Compute T_matrix = D × X and save to disk.  // faster writer + full timing"""
    output_dir = os.path.abspath(output_dir)
    files = os.listdir(output_dir)

    d_files = sorted([f for f in files if re.match(r'Group_\d+_\d+_\d+_D_matrix\.csv', f)])
    x_files = sorted([f for f in files if re.match(r'Group_\d+_\d+_\d+_X_matrix\.csv', f)])

    if len(d_files) != len(x_files):
        logger.warning("Warning: D_matrix and X_matrix file counts do not match.")
        return

    for gi, (d_file, x_file) in enumerate(zip(d_files, x_files), 1):

        df_d = pd.read_csv(os.path.join(output_dir, d_file), memory_map=True)

        df_x = pd.read_csv(os.path.join(output_dir, x_file), memory_map=True)

        d_data = df_d.iloc[:, 1:].to_numpy(dtype=np.int64, copy=False)  # (r, m)
        x_data = df_x.iloc[:, 1:].to_numpy(dtype=np.int64, copy=False)  # (r, n)

        t_matrix = d_data.T.dot(x_data)  # (m, n)

        row_names = list(df_d.columns[1:])
        col_names = list(df_x.columns[1:])
        first_column_name = df_x.columns[0]

        out_path = os.path.join(output_dir, d_file.replace("_D_matrix.csv", "_T_matrix.csv"))

        f = open(out_path, "w", newline="", buffering=4*1024*1024)
        header = first_column_name + "," + ",".join(col_names) + "\n"
        f.write(header)

        try:
            for i, rname in enumerate(row_names):
                f.write(rname)
                f.write(",")
                # key:
                t_matrix[i, :].tofile(f, sep=",", format="%d")
                f.write("\n")
        finally:
            f.close()


def load_encoded_gt_from_matrix_dir(out_dir: str):
    matrix_dir = Path(out_dir) / "matrix_results"

    # Collect all T matrices files
    t_matrix_files = sorted(
        [f for f in matrix_dir.glob("*T_matrix.csv")],
        key=lambda f: (int(f.name.split("_")[1]), int(f.name.split("_")[2]))
    )

    sample_names = []
    gt_buffer = []

    # def encode function
    def encode_line(line, skip_first_col=True):
        parts = line.strip().split(",")
        if skip_first_col:
            parts = parts[1:]
        converted = []
        for val in parts:
            try:
                val = int(val)
                if val == 0:
                    converted.append("0/0")
                elif val == 1:
                    converted.append("1/0")
                elif val >= 2:
                    converted.append("1/1")
                else:
                    converted.append("./.")
            except:
                converted.append("./.")
        return converted

    for idx, file in enumerate(t_matrix_files):
        with open(file, "r") as fin:
            for line_num, line in enumerate(fin):
                if idx == 0 and line_num == 0:
                    # Get the first file's header (sample names)
                    parts = line.strip().split(",")
                    sample_names = parts[1:]  # Drop first column
                    continue
                if line_num == 0:
                    continue  # Other files drop first row
                gt_buffer.append(encode_line(line))
    
    rSV_count = len(gt_buffer)

    return rSV_count, sample_names, gt_buffer


def vcf_generate_from_gt(meta_csv_file: str, output_vcf_file: str, sample_names: list, gt_buffer: list):

    vcf_data = []
    with open(meta_csv_file, newline='') as f_meta:
        reader = list(csv.DictReader(f_meta))
        for index, row in enumerate(reader):
            try:
                group_name = row["group_name"]
                group_parts = group_name.split("_")
                chrom = group_parts[1]
                meta_data = ast.literal_eval(row["meta_array"])
                pos = meta_data.get("pos")
                ref = meta_data.get("ref")
                alt = meta_data.get("alt")
                vcf_data.append([chrom, pos, group_name, ref, alt, ".", ".", "TYPE=rSV", "GT"])
            except Exception as e:
                logger.error(f"Error at row {index}: {e}")
                continue

    with open(output_vcf_file, "w", newline='') as f_out:
        f_out.write("##fileformat=VCFv4.2\n")
        f_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        header_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        f_out.write("\t".join(header_cols + sample_names) + "\n")

        for idx, row in enumerate(vcf_data):
            gt_line = gt_buffer[idx]
            f_out.write("\t".join(map(str, row + gt_line)) + "\n")

    logger.info(f"{output_vcf_file} made successfully.")


def build_rsv_vcf(vcf_dir: str,
                  output_dir: str,
                  meta_csv_file: str,
                  output_vcf_file: str,
                  writedown: bool = False):
    """A switch to control whether to write X/T to disk; otherwise directly output VCF from memory."""
    sample_names, gt_buffer = process_vcf_to_x_matrix(vcf_dir, output_dir, writedown=writedown)

    if writedown:
        # Keep the old path: write X → write T → read T → VCF
        compute_t_matrix(output_dir)
        # Here reuse your original loading function (if your function path is rooted at output_dir
        # instead of matrix_results, please adjust accordingly)
        rSV_count, sample_names2, gt_buffer2 = load_encoded_gt_from_matrix_dir(output_dir)
        assert sample_names2 == sample_names
        vcf_generate_from_gt(meta_csv_file, output_vcf_file, sample_names, gt_buffer2)
    else:
        # New path: directly generate VCF using in-memory GT buffer
        vcf_generate_from_gt(meta_csv_file, output_vcf_file, sample_names, gt_buffer)

    logger.info(f"Done. rSV VCF at: {output_vcf_file}")


# =========================
# Sorting & parsing utilities
# =========================

def _dfile_sort_key_true_style(fname: str) -> Tuple[int, int, int, str]:
    """
    Sort D matrix files by (chrom, group_number, pos, fname).
    File names are like: Group_{chrom}_{group_number}_{pos}_D_matrix.csv
    This matches the sorting of save_meta_csv (chrom, group_number, meta.pos).
    """
    m = re.match(r'Group_(\d+)_(\d+)_(\d+)_D_matrix\.csv', fname)
    if not m:
        return (10**9, 10**9, 10**9, fname)
    chrom, group_number, pos = map(int, m.groups())
    return (chrom, group_number, pos, fname)


def _parse_chrom_group_from_group_name(gname: str) -> Optional[Tuple[int, int]]:
    """
    Extract (chrom, group_number) from group_name.
    Examples:
      'Group_1_42' -> (1, 42)
      '1_42'       -> (1, 42)
    If fewer than two integers are found, return None.
    """
    nums = re.findall(r'\d+', gname)
    if len(nums) < 2:
        return None
    return int(nums[0]), int(nums[1])


def _norm_contig(chrom_token: str, contig_names: List[str]) -> str:
    """
    Normalize chromosome IDs in D files to match VCF header naming.
    - If VCF uses '1,2,3': return the token as is (e.g., '1')
    - If VCF uses 'chr1,chr2': attempt to add 'chr' prefix for matching
    You can also extend this mapping for X/Y/MT if needed.
    """
    if chrom_token in contig_names:
        return chrom_token
    alt = f"chr{chrom_token}"
    if alt in contig_names:
        return alt
    return chrom_token


# =========================
# VCF sequential scan utilities
# =========================

def _make_vcf_iter(vcf_path: str):
    """
    Open an unindexed VCF, return (vcf_in, iterator, contig_order).
    Only a single sequential scan is performed, no tabix random access.
    """
    vcf_in = pysam.VariantFile(vcf_path, "r")
    contigs = list(vcf_in.header.contigs.keys())
    contig_order = {c: i for i, c in enumerate(contigs)}
    return vcf_in, iter(vcf_in), contig_order


def _advance_to(v_iter,
                current,
                target_contig: str,
                target_pos: int,
                contig_order: dict):
    """
    Advance current to the first record at (target_contig, target_pos):
      - If current.contig is before the target, keep advancing;
      - Once at the target contig, advance to pos >= target_pos;
      - Return (current, found), where found indicates if current exactly matches target_pos.
    """
    if current is None:
        return None, False

    # Advance to target contig
    while current and contig_order.get(current.contig, -1) < contig_order.get(target_contig, -1):
        current = next(v_iter, None)

    # Within target contig, advance to >= target_pos
    while current and current.contig == target_contig and current.pos < target_pos:
        current = next(v_iter, None)

    if current and current.contig == target_contig and current.pos == target_pos:
        return current, True
    return current, False


def _read_group_gt_sequential(v_iter,
                              current,
                              contig: str,
                              start_pos: int,
                              n_rows: int,
                              sample_names: List[str],
                              contig_order: dict) -> Tuple[np.ndarray, object]:
    """
    In sequential scan mode, advance from current to (contig, start_pos),
    read n_rows consecutive records, and convert into an (n_rows, n_samples) int16 matrix.
    Values: 0/1/2 (genotype encoding) or -999 (missing/invalid).
    Return (gt_raw, next_current).
    """
    current, ok = _advance_to(v_iter, current, contig, start_pos, contig_order)
    if not ok:
        raise RuntimeError(f"Cannot locate group start at {contig}:{start_pos}. "
                           f"Check D/VCF alignment and sorting.")

    n_samples = len(sample_names)
    A = np.full((n_rows, n_samples), -999, dtype=np.int16)

    for i in range(n_rows):
        if current is None:
            raise RuntimeError(f"VCF ended before reading {n_rows} rows from {contig}:{start_pos}.")
        if current.contig != contig:
            raise RuntimeError(f"Contig switched unexpectedly while reading group at {contig}:{start_pos}.")

        smp = current.samples
        for j, s in enumerate(sample_names):
            gt = smp[s].get("GT")
            if gt is None or len(gt) != 2 or any(a is None for a in gt):
                A[i, j] = -999
            else:
                a, b = gt
                if   (a, b) == (0, 0): A[i, j] = 0
                elif (a, b) in ((0, 1), (1, 0)): A[i, j] = 1
                elif (a, b) == (1, 1): A[i, j] = 2
                else: A[i, j] = -1

        current = next(v_iter, None)

    return A, current


# =========================
# 🔥 Main function: streamed, no large buffer, block-wise rSV.vcf writing
# =========================

def build_rsv_vcf_streamed_nobuf(vcf_dir: str,
                                 output_dir: str,
                                 meta_csv_file: str,
                                 output_vcf_file: str,
                                 strict_meta_group_check: bool = True,
                                 SAMPLE_CHUNK: int = 4096) -> Tuple[int, List[str]]:
    """
    Fully streamed: no X/T writing, no df_vcf, no gt_buffer.
    Strictly follow (chrom, group_number, pos) order of D files,
    aligned with save_meta_csv ordering.
    For each rSV: select original rows by the D column, sum along the sample dimension (in blocks),
    and write to VCF on the fly.

    Returns: (rSV_count, sample_names)
    """
    output_dir = os.path.abspath(output_dir)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # 1) D file order (consistent with save_meta_csv)
    d_files = sorted(
        [f for f in os.listdir(output_dir) if re.match(r'Group_(\d+)_(\d+)_(\d+)_D_matrix\.csv', f)],
        key=_dfile_sort_key_true_style
    )
    if not d_files:
        raise RuntimeError("No D_matrix CSV files found in output_dir.")

    # 2) VCF sequential scan (no tabix)
    vcf_path = os.path.join(vcf_dir, "oSV.vcf")
    vcf_in, v_iter, contig_order = _make_vcf_iter(vcf_path)
    contig_names = list(vcf_in.header.contigs.keys())
    sample_names = list(vcf_in.header.samples)
    n_samples = len(sample_names)

    # 3) Open meta and output VCF (write header first)
    with open(meta_csv_file, newline="") as fmeta, open(output_vcf_file, "w") as fout:
        meta_reader = csv.DictReader(fmeta)

        # VCF Header
        fout.write("##fileformat=VCFv4.2\n")
        fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        fout.write("\t".join(header_cols + sample_names) + "\n")

        rsv_total = 0
        current = next(v_iter, None)

        with tqdm(total=len(d_files), desc="Genotyping rSV (streamed, no-buf)", unit="group", mininterval=0.2) as pbar:
            for dname in d_files:
                m = re.match(r'Group_(\d+)_(\d+)_(\d+)_D_matrix\.csv', dname)
                if not m:
                    pbar.update(1)
                    continue
                chrom_s, group_s, pos_s = m.groups()
                chrom = int(chrom_s)
                group_number = int(group_s)
                start_pos = int(pos_s)

                contig = _norm_contig(str(chrom), contig_names)

                # Read D matrix; first column is block key, following are rSV columns
                df_full = pd.read_csv(os.path.join(output_dir, dname))
                d_mat = df_full.iloc[:, 1:].to_numpy(dtype=np.int8, copy=False)  # small integers save memory
                n_rows, m_rsv = d_mat.shape

                # Sequentially read n_rows lines from VCF at the correct position
                gt_raw, current = _read_group_gt_sequential(
                    v_iter=v_iter,
                    current=current,
                    contig=contig,
                    start_pos=start_pos,
                    n_rows=n_rows,
                    sample_names=sample_names,
                    contig_order=contig_order
                )  # (n_rows, n_samples) int16

                # Write each rSV without constructing full t_mat (m_rsv × n_samples)
                for k in range(m_rsv):
                    try:
                        row = next(meta_reader)
                    except StopIteration:
                        raise RuntimeError("Insufficient meta_csv rows; cannot align with D by group and rSV.")

                    if strict_meta_group_check:
                        got_gname = row.get("group_name", "")
                        parsed = _parse_chrom_group_from_group_name(got_gname)
                        if not parsed:
                            raise RuntimeError(f"Bad group_name in meta: {got_gname}")
                        gchrom, ggroup = parsed
                        if not (gchrom == chrom and ggroup == group_number):
                            raise RuntimeError(
                                "Meta group_name mismatch.\n"
                                f"  D: chrom={chrom}, group={group_number}, start_pos={start_pos}\n"
                                f"  meta group_name={got_gname}\n"
                                "Please ensure meta and D are strictly consistent in (chrom, group_number) "
                                "(pos should be sorted within meta)."
                            )

                    meta = ast.literal_eval(row["meta_array"])
                    pos = int(meta.get("pos"))
                    ref = str(meta.get("ref"))
                    alt = str(meta.get("alt"))
                    group_name = row["group_name"]
                    info = "TYPE=rSV"

                    # Original rows participating in this rSV (usually sparse)
                    rows_idx = np.flatnonzero(d_mat[:, k])

                    # Write the fixed 9 columns
                    fout.write("\t".join([contig, str(pos), group_name, ref, alt, ".", ".", info, "GT"]))

                    if rows_idx.size == 0:
                        # No contributing rows: all 0/0 (equivalent to empty sum)
                        fout.write("\t" + "\t".join(("0/0",) * n_samples) + "\n")
                        rsv_total += 1
                        continue

                    # Process samples in chunks to avoid building oversized strings
                    j0 = 0
                    while j0 < n_samples:
                        j1 = min(j0 + SAMPLE_CHUNK, n_samples)
                        # Sum genotypes: treat -999 as invalid (not 0/1/2), encode as missing later
                        s = gt_raw[rows_idx, j0:j1].astype(np.int32, copy=False).sum(axis=0)

                        enc_chunk = []
                        for val in s:
                            if val == 0:
                                enc_chunk.append("0/0")
                            elif val == 1:
                                enc_chunk.append("1/0")
                            elif val >= 2:
                                enc_chunk.append("1/1")
                            else:
                                enc_chunk.append("./.")
                        fout.write("\t" + "\t".join(enc_chunk))
                        j0 = j1

                    fout.write("\n")
                    rsv_total += 1

                pbar.update(1)

    logger.info(f"rSV.vcf written: {output_vcf_file}")
    return rsv_total, sample_names
