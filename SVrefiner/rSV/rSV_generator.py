import os
import re
import click
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from pathlib import Path
from collections import OrderedDict
from ..utils.logging_utils import get_logger, log_tqdm_summary
from concurrent.futures import ThreadPoolExecutor, as_completed


logger = get_logger(__name__)


# === Encoding Map ===
base_map = {"-": 0, "A": 1, "a": 1, "T": 2, "t": 2, "C": 3, "c": 3, "G": 4, "g": 4, "N": 5, "n": 5}
reverse_map = {0: "-", 1: "A", 2: "T", 3: "C", 4: "G", 5: "N"}


def encode_sequence(seq):
    return np.array([base_map.get(base.upper(), 5) for base in seq], dtype=int)


def num_to_base(n):
    return reverse_map.get(n, "-")


# === FASTA Input Utilities ===

def parse_filename_metadata(filename):
    pattern = r"Group_(\d+)_(\d+)_(\d+)_aligned\.fasta"
    match = re.search(pattern, filename)
    if match:
        return {
            "chrom": int(match.group(1)),
            "group_num": int(match.group(2)),
            "pos": int(match.group(3))
        }
    return None


def list_aligned_fasta_files(folder):
    aligned_files = []
    if not os.path.exists(folder):
        click.echo(f"[Warning] Alignment folder not found: {folder}")
        return []

    for fname in os.listdir(folder):
        if fname.endswith("_aligned.fasta"):
            full_path = os.path.join(folder, fname)
            meta = parse_filename_metadata(fname)
            if meta:
                aligned_files.append({"path": full_path, "meta": meta})

    click.echo(f"Detected {len(aligned_files)} aligned FASTA files.")
    return aligned_files


def preload_fasta_files(file_infos):
    in_memory_data = []
    for file_info in file_infos:
        records = list(SeqIO.parse(file_info["path"], "fasta"))
        info_copy = file_info.copy()
        info_copy["records"] = records
        in_memory_data.append(info_copy)
    click.echo(f"Loaded {len(in_memory_data)} FASTA files into memory.")
    return in_memory_data


# === Block Construction Utilities ===

def is_rank1_block(mat):
    if np.all(mat == 0):
        return False
    return np.linalg.matrix_rank(mat) == 1


def is_rank1_pair(col1, col2):
    return np.linalg.matrix_rank(np.stack([col1, col2], axis=1)) == 1


def find_blocks_with_rank(matrix, window=10):
    blocks = []
    i = 0
    while i < matrix.shape[1]:
        if np.all(matrix[:, i] == 0):
            blocks.append((i, i))
            i += 1
            continue

        end = min(i + window, matrix.shape[1])
        submat = matrix[:, i:end]
        if is_rank1_block(submat):
            blocks.append((i, end - 1))
            i = end
        else:
            start = i
            i += 1
            while i < matrix.shape[1]:
                if np.all(matrix[:, i] == 0):
                    blocks.append((start, i - 1))
                    blocks.append((i, i))
                    i += 1
                    break
                if not is_rank1_pair(matrix[:, i - 1], matrix[:, i]):
                    blocks.append((start, i - 1))
                    start = i
                i += 1
            else:
                if start < matrix.shape[1]:
                    blocks.append((start, matrix.shape[1] - 1))
    return blocks


def find_blocks_mask_based(matrix):
    """
    Identify column blocks where adjacent columns share the same nonzero mask.
    
    A block is defined as a consecutive set of columns where each pair of 
    adjacent columns has exactly the same nonzero positions.
    
    Args:
        matrix (np.ndarray): Input matrix of shape (m, n), where columns represent vectors.
    
    Returns:
        list of tuples: Each tuple (start, end) represents the column index range of a block.
    """
    # Generate a boolean mask where True indicates nonzero entries
    mask = matrix != 0  # shape: (m, n)

    # Check if adjacent columns share the same nonzero mask
    is_rank1 = np.all(mask[:, :-1] == mask[:, 1:], axis=0)  # shape: (n - 1,)

    # Identify blocks based on the pairwise mask comparison
    blocks = []
    start = 0
    for i, related in enumerate(is_rank1):
        if not related:
            blocks.append((start, i))
            start = i + 1
    blocks.append((start, matrix.shape[1] - 1))  # Add the final block

    return blocks


def smart_split_by_row_patterns(mat):
    unique_patterns, inverse_indices = np.unique(mat, axis=0, return_inverse=True)
    sub_matrices = []
    for idx in range(unique_patterns.shape[0]):
        mask = (inverse_indices == idx)
        sub_mat = np.zeros_like(mat)
        sub_mat[mask] = mat[mask]
        if np.any(sub_mat != 0):
            sub_matrices.append(sub_mat)
    return sub_matrices


def extract_blocks_with_split(matrix, blocks, start_index_shift=True):
    expanded, new_blocks = [], []
    for (start, end) in blocks:
        mat = matrix[:, start:end + 1]
        if np.unique(mat, axis=0).shape[0] <= 2:
            expanded.append(mat)
            new_blocks.append((start, end))
        else:
            sub_mats = smart_split_by_row_patterns(mat)
            for i, submat in enumerate(sub_mats):
                new_start = start + i * mat.shape[1] if start_index_shift else 0
                new_end = new_start + mat.shape[1] - 1
                expanded.append(submat)
                new_blocks.append((new_start, new_end))
    return expanded, new_blocks


# === Matrix and Metadata Construction ===

def convert(row):
    return ''.join([num_to_base(x) for x in row if x != 0])


def build_D_and_meta(expanded_mats):
    D_cols, metas = [], []
    ref_prefix, pos, pos_buffer = "", 0, 0

    for i, mat in enumerate(expanded_mats):
        first_col = mat[:, 0]
        last_col = mat[:, -1]

        if i == 0:
            ref_prefix = num_to_base(last_col[0])
            continue

        d_col = (first_col != 0).astype(int) if first_col[0] == 0 else (first_col == 0).astype(int)
        D_cols.append(d_col[:, None])

        row_idx = np.where(d_col == 1)[0]
        if len(row_idx) == 0:
            continue

        ref = ref_prefix + convert(mat[0])
        alt = ref_prefix + convert(mat[row_idx[0]])

        pos_buffer = pos
        
        if last_col[0] != 0:
            ref_prefix = num_to_base(last_col[0])
            
            pos += mat.shape[1]

        metas.append({"pos": pos_buffer, "ref": ref, "alt": alt})

    D = np.hstack(D_cols)[1:, :]
    return D, metas


def merge_rsv_by_pos_and_mask(D: np.ndarray, metas: list):
    """
    Only used for the has_insertion group:
    - Within the same group, merge columns where meta['pos'] is identical,
      and further merge if the column masks from D are identical.
    - Merge rule: ref takes the first column of the group;
      alt = alt1 + alt2[1:] + alt3[1:] + ...
    - The D matrix only keeps the first column of each identical mask group
      (since the masks are the same).
    - Returns: D_new (ndarray), metas_new (list[dict]), with columns ordered
      by ascending pos and preserving the original order within groups.
    """
    if D.size == 0 or len(metas) == 0:
        return D, metas

    # 1) Group column indices by pos (preserve stable order)
    pos_to_cols = OrderedDict()
    for j, m in enumerate(metas):
        p = int(m["pos"])
        pos_to_cols.setdefault(p, []).append(j)

    keep_cols = []
    metas_new = []

    # 2) For each pos group, further group by mask
    for p in sorted(pos_to_cols.keys()):
        cols = pos_to_cols[p]  # columns with the same pos
        mask_groups = OrderedDict()  # mask(tuple) -> list of column indices
        for j in cols:
            mask = tuple(int(x) for x in D[:, j].tolist())
            mask_groups.setdefault(mask, []).append(j)

        # 3) Merge each mask group: keep the first column, concatenate alt[1:]
        for mask, js in mask_groups.items():
            js_sorted = sorted(js)            # keep original column order
            j0 = js_sorted[0]                 # representative column of this mask group
            merged_ref = metas[j0]["ref"]
            merged_alt = metas[j0]["alt"]
            for j in js_sorted[1:]:
                aj = metas[j]["alt"]
                merged_alt += (aj[1:] if len(aj) > 0 else "")
            keep_cols.append(j0)
            metas_new.append({"pos": p, "ref": merged_ref, "alt": merged_alt})

    # 4) Reconstruct new D (columns follow metas_new order)
    D_new = D[:, keep_cols] if len(keep_cols) > 0 else D

    return D_new, metas_new


def process_fasta_in_memory(file_info, has_insertion_dict):
    records = file_info["records"]
    pos_prefix = file_info["meta"]["pos"]
    fasta_path = file_info["path"]
    group_name = re.sub(r'_aligned\.fasta$', '', os.path.basename(fasta_path))

    seqs = [str(rec.seq) for rec in records]
    matrix = np.vstack([encode_sequence(seq) for seq in seqs])

    blocks = find_blocks_with_rank(matrix) if has_insertion_dict[group_name] == True else find_blocks_mask_based(matrix)
    expanded_mats, _ = extract_blocks_with_split(matrix, blocks)

    D, metas = build_D_and_meta(expanded_mats)

    for m in metas:
        m["pos"] += pos_prefix

    if has_insertion_dict.get(group_name, False):
        D, metas = merge_rsv_by_pos_and_mask(D, metas)

    return os.path.basename(fasta_path), metas, D


# === Output ===

def save_D_matrix(D, seq_ids, sv_ids_full, out_path):
    sv_ids = [x.split("_")[-1] for x in sv_ids_full]
    df = pd.DataFrame(D, columns=sv_ids)
    df.insert(0, "sv_id", seq_ids) 
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, index=False)


def save_meta_csv(meta_list, out_path):
    meta_tuples = [(item["group_name"], eval(item["meta_array"])) for item in meta_list]
    meta_tuples.sort(key=lambda x: (int(x[0].split('_')[1]), int(x[0].split('_')[2]), x[1]["pos"]))
    sorted_meta = [{"group_name": g, "meta_array": str(m)} for g, m in meta_tuples]
    df = pd.DataFrame(sorted_meta)
    df.to_csv(out_path, index=False)


# === Main Pipeline ===

def rSV_generator(out: str, has_insertion_dict: dict, threads: int = 16):
    matrix_dir = Path(out) / "matrix_results"
    matrix_dir.mkdir(parents=True, exist_ok=True)

    alignments_dir = Path(out) / "alignment_results"
    all_files = list_aligned_fasta_files(alignments_dir)
    in_memory_data = preload_fasta_files(all_files)

    results = []
    failure_count = 0

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_fasta_in_memory, item, has_insertion_dict): item for item in in_memory_data}
        pbar = tqdm(as_completed(futures), total=len(futures), desc="rSV generating", unit="groups")
        for future in pbar:
            try:
                fname, metas, D = future.result()
                meta_info = futures[future]["meta"]
                results.append({
                    "source": fname,
                    "chrom": meta_info["chrom"],
                    "group": meta_info["group_num"],
                    "original_pos": meta_info["pos"],
                    "meta": metas,
                    "D": D
                })
            except Exception as e:
                logger.error(f"Failed to process file: {futures[future]['path']} — {repr(e)}")
                failure_count += 1
        log_tqdm_summary(pbar, logger)

    meta_records = []
    for item in results:
        fname = item["source"]
        metas = item["meta"]
        D = item["D"]
        group_name = fname.replace("_aligned.fasta", "")
        sv_ids_full = [f"{group_name}_rSV{i+1}" for i in range(D.shape[1])]
        seq_ids = [f"seq{i+1}" for i in range(D.shape[0])]
        save_D_matrix(D, seq_ids, sv_ids_full, matrix_dir / f"{group_name}_D_matrix.csv")
        for i, meta in enumerate(metas, start=1):
            meta_records.append({
                "group_name": f"{group_name}_rSV{i}",
                "meta_array": str(meta)
            })

    try:
        save_meta_csv(meta_records, Path(out) / "rSV_meta.csv")
    except Exception as e:
        logger.error(f"Failed to save rSV_meta.csv: {repr(e)}")
