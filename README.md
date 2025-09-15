
# SVrefiner

**SVrefiner** is a Python-based tool for processing structural variants (SVs) and generating refined SVs (rSVs) in VCF v4.2 format.
It integrates sequence alignment, window scanning, and merging workflows to resolve overlapping SVs, with a focus on `DEL` (deletions) and `INS` (insertions).
Complex SVs (ref* and alt* differ at the first base, or both ref and alt are longer than 1 bp) are not processed in the current version to improve efficiency, as their occurrence in real datasets is relatively small (e.g., only **2.2%** in the tomato pangenome; Zhou *et al.*, *Nature*, 2022, 606: 527–534).
> ref: reference allele, alt: alternative allele
---

## Quick Start

### Installation (Conda recommended)

```bash
# Create and activate environment
conda create -n SVrefiner python=3.8
conda activate SVrefiner

# Clone and install SVrefiner
git clone https://github.com/Lostmet/SVrefiner.git
cd SVrefiner
pip install .

# Install MAFFT
conda install conda-forge::mafft
```

---

## Input

**Required input files:**

* **VCF** file with index (`.vcf.gz` and `.vcf.gz.tbi` or `.vcf.gz.csi`)
* **FASTA** reference genome (`.fa` or `.fasta`) including the chromosomes covered by SVs

---

## Main Output

| File/Directory          | Description                                   |
| ----------------------- | --------------------------------------------- |
| `alignment_results/`    | Alignment results for each overlapping group  |
| `matrix_results/`       | Output matrices for each overlapping group    |
| `rSV.vcf`               | Refined SVs (final VCF)                       |
| `*.log`                 | Execution logs                                |

---

### Notes

- **Chromosome identifiers in VCF and FASTA must be numeric.**
  - Examples: `1`, `2`, `3` (autosomes), `23` (X), `24` (Y)
  - Identifiers such as `chr1`, `chrX` are not supported.
  
- **SVs must be normalized and left-aligned in standard VCF format.**
  - Only **biallelic SVs** are supported; multiallelic records are not accepted.
  - If your VCF is not normalized, please preprocess it using external tools (e.g., `bcftools norm`) before running SVrefiner.
  
- It is recommended to perform quality control (QC) on the SVs prior to running SVrefiner.
  - Low-quality SVs, especially those with high missingness (low call rate), can significantly compromise the quality of the resulting rSVs, please consider filtering such variants before proceeding.


---

## VCF Example

> The VCF file must be compressed (`.vcf.gz`) and indexed (`.tbi` or `.csi`).

```vcf
##fileformat=VCFv4.2
##source=YourTool
#CHROM  POS  ID    REF     ALT     QUAL    FILTER  INFO   FORMAT Sample1  Sample2  Sample3  Sample4
1       1    sv1   ACTA    A       50      PASS    .       GT      1/1      1/0      0/0      ./.
1       5    sv2   G       GAAC    99      PASS    .       GT      0/0      1/0      0/0      0/0
1       6    sv3   GCTAG   <INV>   98      PASS    .       GT      ./.      0/0      1/1      1/1
```

---

## Commands

| Command       | Function                                                                |
| ------------- | ----------------------------------------------------------------------- |
| `process-vcf` | Process VCF to identify overlapping SVs; output `nSV.vcf` and `oSV.vcf` |
| `align`       | Perform sequence alignment for overlapping SVs                          |
| `make-rsv`    | Define refined SVs from aligned sequences and output `rSV.vcf`          |
| `run-all`     | Execute the complete pipeline                                           |

## Options
| Option           | Description                                                                         |
| ---------------- | ----------------------------------------------------------------------------------- |
| `--vcf`          | **Required.** Input VCF file (`.vcf.gz`; index `.csi` or `.tbi` required).          |
| `--ref`          | **Required.** Input reference FASTA file (`.fasta` or `.fa`).                       |
| `--out`          | **Required.** Output directory.                                                     |
| `--threads`      | *Optional.* Number of threads to use (default: `10`).                               |
| `--write-matrix` | *Optional.* `YES`/`NO` (default: `NO`). Whether to output **X** and **T** matrices. |


### Example Workflow

#### Step 1: Process VCF

```bash
SVrefiner process-vcf --vcf test.vcf.gz --ref test.fasta --out test --threads 10
```

This separates overlapping (`oSV.vcf`) and non-overlapping (`nSV.vcf`) SVs.

#### Step 2: Align

```bash
SVrefiner align --vcf test.vcf.gz --ref test.fasta --out test --threads 10
```

#### Step 3: Define rSV

```bash
SVrefiner make-rsv --vcf test.vcf.gz --ref test.fasta --out test --threads 10
```

#### Complete Pipeline

```bash
SVrefiner run-all --vcf test.vcf.gz --ref test.fasta --out test --threads 10
```

---

## Requirements

* **Python ≥ 3.8**
* **MAFFT ≥ v7.526**
* Python libraries:

  * `pandas`
  * `numpy`
  * `biopython`
  * `click`
  * `tqdm`

---

## Conceptual Overview

The pipeline resolves overlapping SVs into more precise rSVs through alignment and window-based clustering (Figure 1).

<p align="center">
<img src="https://github.com/user-attachments/assets/69c43992-af28-4669-a35e-198771986241" width="800">
</p>
<p align="center"><b>Figure 1.</b> Schematic of SVrefiner workflow and outputs.</p>

For detailed algorithmic description, see: [Algorithmic Logic of rSV Software](https://github.com/Lostmet/Algorithmic_Logic_of_rSV_Software).

---

## Output Structure

### Main Directory

| File/Directory | Description         |
| -------------- | ------------------- |
| `rSV.vcf`      | Final refined SVs   |
| `rSV_meta.csv` | rSV metadata        |
| `nSV.vcf`      | Non-overlapping SVs |
| `oSV.vcf`      | Overlapping SVs     |
| `*.log`        | Runtime logs        |

Example log:

```
Command:  SVrefiner run-all
  --ref   test.fasta 
  --vcf   test.vcf.gz 
  --out   test 
Start time:                   2025-07-07 13:24:57
End time:                     2025-07-07 13:24:59
Total runtime:                0:00:01

Total variants:               100
INV count:                    1
nSV count:                    72
Excluded SNP count:           0
Overlapping SVs:              28
Overlapping SVs percentage:   28.00%
Total variant groups:         9
Final rSV count:              45
```

### `alignment_results/`

| File                          | Description                |
| ----------------------------- | -------------------------- |
| `Group_input_origin.fasta`    | Original group sequences   |
| `Group_input_spliced.fasta`   | Sliced insertion sequences |
| `Group_aligned_spliced.fasta` | Aligned slices             |
| `Group_aligned.fasta`         | Aligned sequences          |

### `matrix_results/`

| File                 | Description                         |
| -------------------- | ----------------------------------- |
| `Group_D_matrix.csv` | D matrix (rSV–SV relationships)     |
| `Group_X_matrix.csv` | X matrix (SV–sample relationships)  |
| `Group_T_matrix.csv` | T matrix (rSV–sample relationships) |

### `alignment_error_logs/`

Contains MAFFT error logs, if any.

---

## Support

If you encounter issues, please open an issue on GitHub or contact the authors at:

* [fenglostmet@tju.edu.cn](mailto:fenglostmet@tju.edu.cn)
* [xia\_xiaoxuan@outlook.com](mailto:xia_xiaoxuan@outlook.com)

---
## Citation

If you use this software, please cite:

Xia X., Wu J., Gao Z. *et al.*, Modeling structural variations sequencing information to address missing heritability and enhance risk prediction (2025) *bioRxiv*. doi:[10.1101/2025.08.07.669060](https://doi.org/10.1101/2025.08.07.669060)




