import os
import time
import click
import shutil
import logging
import datetime
from pathlib import Path
from .config import Config
from .utils.logging_utils import *
from .genotype.genotype_mapper import *
from .genotype.matrix_generator import *
from .rSV.rSV_generator import rSV_generator
from .alignment.alignment import run_alignments
from .utils.file_utils import check_vcf_vs_fasta
from .variant_processing.fasta_generator import generate_fasta_sequences
from .variant_processing.vcf_parser import process_variants, filter_vcf, make_oSV


logger = get_logger(__name__)


def check_tools(*tools):
    """
    Check if required tools are available in the user's PATH.
    Raise an error if any tool is missing.
    """
    missing_tools = [tool for tool in tools if shutil.which(tool) is None]
    if missing_tools:
        raise RuntimeError(
            f"The following tools are missing from your environment: {', '.join(missing_tools)}. "
            f"Please ensure they are installed and accessible from your PATH."
        )

@click.group()
@click.version_option("0.1.0", prog_name="SVrefiner", message="""
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
SVrefiner - A Python tool for generating refined Structured Variants (rSVs)  
Version 0.1.0 - Contact: fenglostmet@tju.edu.cn
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
""")
def cli():
    """A Python tool for generating refined Structured Variants (rSVs)."""
    pass


@cli.command("align")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory for processed variants and FASTA files')
@click.option('--threads', default=1, type=int, help='Number of threads')
def align_time(vcf: str, ref: str, out: str, threads: int):
    """Do alignment of overlapping SVs"""
    setup_logging(Path(out) / "align.log")
    start_time = time.time()
    log_step("Step 1 Processing VCF and generating FASTA")
    config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads=threads)
    grouped_variants_list, _, _, _, \
    _, _, _, _ = process_variants(config)

    grouped_variants_dict = {}
    for group in grouped_variants_list:
        grouped_variants_dict.setdefault(group.chrom, []).append(group)
    total_groups = sum(len(groups) for groups in grouped_variants_dict.values())
    logger.info(f"Overlapping variants grouped: {total_groups:,} groups")

    fasta_path, has_insertion_dict, poly_ins_list = generate_fasta_sequences(config, grouped_variants_dict, total_groups)
    logger.info(f"FASTA file created: {fasta_path}")

    log_step("Step 2 Running alignments")
    alignments_config = Config(grouped_variants_file=fasta_path, ref_fasta=ref, output_dir=out, threads=threads)
    run_alignments(alignments_config, fasta_path, has_insertion_dict, poly_ins_list)
    logger.info("Alignments completed. Results saved in 'alignment_results' directory.")    
    end_time = time.time()
    total_time = end_time - start_time
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    log_step(f"Total runtime: {int(hours)}:{int(minutes):02d}:{int(seconds):02d}")

@cli.command("process-vcf")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory for processed variants and FASTA files')
@click.option('--threads', default=8, type=int, help='Number of threads')

def process_vcf(vcf: str, ref: str, out: str, threads: int):
    """Process overlapping variants, group them, and generate nSVs, oSVs. 
    This step produces 'variants_pre_aligned.fasta', 'oSV,vcf' and 'nSV.vcf'."""
    setup_logging(Path(out) / "process_vcf.log")
    try:
        logger.info(f"{'Command:':<5}{get_clean_command()}")
        start_time = time.time()
        logger.info(f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        log_thread_info(threads)

        log_step("Processing VCF and generating FASTA")
        config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads=threads)
        grouped_variants_list, single_sv_count, \
        multi_bp, percentage_sv_overlapped, \
        single_group, inv_count, variant_count, snp_groups = process_variants(config)

        est_time = multi_bp / 3000
        hours, remainder = divmod(est_time, 3600)  
        minutes, _ = divmod(remainder, 60)  
        click.echo(f"Estimated runtime: {int(hours)}:{int(minutes):02d}:00")
        click.echo(click.style("Note: Actual runtime depends on CPU performance.", fg="yellow"))

        click.echo("Processing non-overlapping SVs (nSVs) and generating nSV.vcf...")
        _ = filter_vcf(config, single_group, snp_groups)

        click.echo("Processing overlapping SVs (oSVs) and generating oSV.vcf...")
        make_oSV(config, grouped_variants_list) # make oSV.vcf
        
        stats = {
                    "Total variants": f"{variant_count:,}",
                    "INV count": f"{inv_count:,}",
                    "nSV count": f"{single_sv_count:,}",
                    "Excluded SNP count": f"{len(snp_groups)}",
                    "Overlapping SVs": f"{variant_count - single_sv_count:,}",
                    "Overlapping SVs percentage": f"{percentage_sv_overlapped:.2f}%",
                    "Total variant groups": f"{len(grouped_variants_list):,}",
                }
  
        end_time = time.time()
        runtime = end_time - start_time
        log_step("Summary")
        log_summary_block(
            cmd=get_clean_command(),
            start=start_time,
            duration=runtime,
            stats=stats)
        log_all_warnings_and_errors()

    except Exception as e:
        logger.info(f"Error in process-vcf: {str(e)}", err=True)
        logging.error(f"Error in process-vcf: {str(e)}")
        raise click.Abort()


@cli.command("run-all")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory')
@click.option('--threads', default=10, type=int, help='Number of threads')

def run_all(vcf: str, ref: str, out: str, threads: int):
    """Execute the complete variant processing pipeline."""
    setup_logging(Path(out) / "run-all.log")
    check_tools("mafft")
    check_vcf_vs_fasta(vcf, ref)

    try:
        logger.info(f"{'Command:':<5}{get_clean_command()}")
        start_time = time.time()
        logger.info(f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        log_thread_info(threads)

        log_step("Step 1 Processing VCF and generating FASTA")
        config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads=threads)
        grouped_variants_list, single_sv_count, \
        multi_bp, percentage_sv_overlapped, \
        single_group, inv_count, variant_count, snp_groups = process_variants(config)

        est_time = multi_bp / 3000
        hours, remainder = divmod(est_time, 3600)  
        minutes, _ = divmod(remainder, 60)  
        click.echo(f"Estimated runtime: {int(hours)}:{int(minutes):02d}:00")
        click.echo(click.style("Note: Actual runtime depends on CPU performance.", fg="yellow"))

        click.echo("Processing non-overlapping SVs (nSVs) and generating nSV.vcf...")
        nSV_name = filter_vcf(config, single_group, snp_groups)

        click.echo("Processing overlapping SVs (oSVs) and generating oSV.vcf...")
        make_oSV(config, grouped_variants_list) # make oSV.vcf
        
        grouped_variants_dict = {}
        for group in grouped_variants_list:
            grouped_variants_dict.setdefault(group.chrom, []).append(group)
        total_groups = sum(len(groups) for groups in grouped_variants_dict.values())
        logger.info(f"Overlapping variants grouped: {total_groups:,} groups")

        fasta_path, has_insertion_dict, poly_ins_list = generate_fasta_sequences(config, grouped_variants_dict, total_groups)
        logger.info(f"FASTA file created: {fasta_path}")

        log_step("Step 2 Running alignments")
        alignments_config = Config(grouped_variants_file=fasta_path, ref_fasta=ref, output_dir=out, threads=threads)
        run_alignments(alignments_config, fasta_path, has_insertion_dict, poly_ins_list)
        logger.info("Alignments completed. Results saved in 'alignment_results' directory.")

        log_step("Step 3 Generating rSVs")
        matrix_dir = Path(out) / "matrix_results"
        rSV_generator(out, has_insertion_dict, threads)

        click.echo("Creating X matrices...")
        sample_names = process_vcf_to_x_matrix(out, matrix_dir)

        click.echo("Creating T matrices...")
        compute_t_matrix(matrix_dir)

        logger.info("Matrices successfully generated in 'matrix_results'.")

        rSV_meta_csv = os.path.join(out, "rSV_meta.csv")
        gt_matrix = os.path.join(out, "GT_matrix.csv")
        rSV_vcf = os.path.join(out, "rSV.vcf")

        click.echo("Generating GT matrix and rSV.vcf...")
        rSV_count, sample_names, gt_buffer = load_encoded_gt_from_matrix_dir(out)
        logging.info(f"rSV count: {rSV_count:,}")
        vcf_generate_from_gt(rSV_meta_csv, rSV_vcf, sample_names, gt_buffer)

        logging.info(f"{rSV_vcf} successfully generated!")
        click.echo("Pipeline execution completed!")

        stats = {
            "Total variants": f"{variant_count:,}",
            "INV count": f"{inv_count:,}",
            "nSV count": f"{single_sv_count:,}",
            "Excluded SNP count": f"{len(snp_groups)}",
            "Overlapping SVs": f"{variant_count - single_sv_count:,}",
            "Overlapping SVs percentage": f"{percentage_sv_overlapped:.2f}%",
            "Total variant groups": f"{total_groups:,}",
            "Final rSV count": f"{rSV_count:,}"
        }
  
        end_time = time.time()
        runtime = end_time - start_time
        log_step("Summary")
        log_summary_block(
            cmd=get_clean_command(),
            start=start_time,
            duration=runtime,
            stats=stats)
        log_all_warnings_and_errors()
    except RuntimeError as e:
        logger.warning(f"Error: {str(e)}", err=True)
        raise click.Abort()
    except Exception as e:
        logger.warning(f"Critical error in make-rSV: {str(e)}", err=True)
        raise click.Abort()


@cli.command("make-rsv")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--out', required=True, help='Output directory')
@click.option('--ref', default=None, help='Reference FASTA file')
@click.option('--threads', default=10, type=int, help='Number of threads')

def make_rSV(vcf: str, ref, out: str, threads: int):
    '''If you already have alignment results and oSV.vcf, you can make rSVs by this tag. 
    Ensure the output path remains unchanged.'''
    setup_logging(Path(out) / "make-rsv.log")
    start_time = time.time()
    log_step("Step 1 Extracting SV information")
    config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads=threads)

    grouped_variants_list, single_sv_count, \
    _, percentage_sv_overlapped, \
    _, inv_count, variant_count, snp_groups = process_variants(config)

    log_step("Step 2 Converting to VCF format and generating matrices")

    grouped_variants_dict = {}
    for group in grouped_variants_list:
        grouped_variants_dict.setdefault(group.chrom, []).append(group)
    total_groups = sum(len(groups) for groups in grouped_variants_dict.values())
    logger.info(f"Overlapping variants grouped: {total_groups:,} groups")

    fasta_path, has_insertion_dict, poly_ins_list = generate_fasta_sequences(config, grouped_variants_dict, total_groups)

    # Create a directory for storing matrices
    matrix_dir = Path(out) / "matrix_results"
    
    sample_names = sample_name_contract(vcf)
    
    click.echo("1. Generating rSV_meta.csv with metadata (positions, reference, and alternate alleles)...")
    rSV_generator(out, has_insertion_dict, threads)

    click.echo("2. Creating matrices...")
    sample_names = process_vcf_to_x_matrix(out, matrix_dir)
    compute_t_matrix(matrix_dir)
    logger.info("D, X, and T matrices successfully generated in 'matrix_results' directory.")
    
    rSV_meta_csv = os.path.join(out, "rSV_meta.csv")
    gt_matrix = os.path.join(out, "GT_matrix.csv")
    rSV_vcf = os.path.join(out, "rSV.vcf")
    
    click.echo("3. Generating GT matrix and rSV.vcf...")
    rSV_count, sample_names, gt_buffer = load_encoded_gt_from_matrix_dir(out)
    logging.info(f"rSV count: {rSV_count:,}")
    vcf_generate_from_gt(rSV_meta_csv, rSV_vcf, sample_names, gt_buffer)

    logger.info(f"{rSV_vcf} successfully generated!")
    click.echo("Pipeline execution completed!")
    click.echo("All steps completed successfully!")

    stats = {
        "Total variants": f"{variant_count:,}",
        "INV count": f"{inv_count:,}",
        "nSV count": f"{single_sv_count:,}",
        "Excluded SNP count": f"{len(snp_groups)}",
        "Overlapping SVs": f"{variant_count - single_sv_count:,}",
        "Overlapping SVs percentage": f"{percentage_sv_overlapped:.2f}%",
        "Total variant groups": f"{total_groups:,}",
        "Final rSV count": f"{rSV_count:,}"
    }

    end_time = time.time()
    runtime = end_time - start_time
    log_step("Summary")
    log_summary_block(
        cmd=get_clean_command(),
        start=start_time,
        duration=runtime,
        stats=stats)
    log_all_warnings_and_errors()
