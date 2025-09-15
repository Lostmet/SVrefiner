import pysam
import logging
import click
import re

def check_vcf_vs_fasta(vcf_path: str, ref_path: str):
    logger = logging.getLogger()

    # STEP 1: Contract chrome ID of VCF file
    try:
        with pysam.VariantFile(vcf_path) as vcf:
            vcf_chroms = set(rec.chrom for rec in vcf.fetch())
    except Exception as e:
        logger.error(f"Failed to open or parse VCF file: {e}")
        raise click.Abort()

    if not vcf_chroms:
        logger.warning("No variants found in VCF. Skipping chromosome consistency check.")
        return

    # STEP 2: Detect "chr"
    if any(chrom.startswith("chr") for chrom in vcf_chroms):
        logger.warning("Detected 'chr' prefix in chromosome names. Please encode them as integers.")
        raise click.Abort()

    # STEP 3: Check chromosome ID
    non_integer_chroms = sorted([c for c in vcf_chroms if not re.fullmatch(r"\d+", c)])
    if non_integer_chroms:
        logger.error("Detected non-integer chromosome names in VCF:")
        for chrom in non_integer_chroms:
            logger.error(f"  - {chrom}")
        logger.error("Please encode all chromosome names as integers.")
        raise click.Abort()

    # STEP 4: Contract chromosome ID of FASTA file
    try:
        ref = pysam.FastaFile(ref_path)
        ref_chroms = set(ref.references)
    except Exception as e:
        logger.error(f"Failed to open reference FASTA file: {e}")
        raise click.Abort()

    # STEP 5: Compare with VCF to make sure that matches
    missing = sorted(chrom for chrom in vcf_chroms if chrom not in ref_chroms)
    if missing:
        logger.error("The following chromosomes are present in the VCF but missing from the reference FASTA:")
        for chrom in missing:
            logger.error(f"  - {chrom} (expected '> {chrom}' in FASTA)")

        logger.error("Chromosomes available in the reference FASTA:")
        for ref_name in sorted(ref_chroms):
            logger.error(f"  > {ref_name}")

        logger.error("Aborting due to inconsistent chromosome naming between VCF and reference FASTA.")
        raise click.Abort()
