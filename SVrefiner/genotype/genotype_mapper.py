import ast
import csv
import sys
from pathlib import Path
from ..utils.logging_utils import get_logger


logger = get_logger(__name__)
csv.field_size_limit(sys.maxsize)


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
