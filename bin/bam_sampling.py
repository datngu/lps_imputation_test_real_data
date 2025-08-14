#!/usr/bin/env python3
import argparse
import subprocess
import random
import pysam

# Reproducibility
SEED = 2024
random.seed(SEED)

# Genome size reference (GRCh38/GRCh37 ~ 3.3 Gbp)
GENOME_SIZE = 3.3e9
SD_COV = 0.1              # Standard deviation of coverage
MIN_COV = 0.1             # Minimum coverage
READ_LENGTH = 150         # Illumina 150 bp
DUPLICATED_RATE = 0.09    # Estimated duplication rate

def depth2count(target_depth):
    """
    Convert target mean depth into expected read count considering duplicates and variation.
    """
    depth = (1 - DUPLICATED_RATE) * target_depth
    min_cut = max(MIN_COV, depth - 1.28 * SD_COV)
    max_cut = depth + 1.28 * SD_COV

    ran_depth = random.gauss(depth, SD_COV)
    while not (min_cut <= ran_depth <= max_cut):
        ran_depth = random.gauss(depth, SD_COV)

    # Number of reads needed
    n_read = (ran_depth * GENOME_SIZE) // READ_LENGTH
    return int(n_read)

def count_total_reads(bam_path):
    """
    Count total mapped reads using BAM index statistics.
    """
    subprocess.run(["samtools", "index", bam_path], check=True)
    bam = pysam.AlignmentFile(bam_path, "rb")
    return sum(stat[3] for stat in bam.get_index_statistics())



def run_picard_downsample(input_bam, output_bam, proportion):
    """
    Run Picard DownsampleSam via Java command.
    """
    cmd = [
        "java", "-Xmx32g", "-jar", "/biotools/picard.jar",
        "DownsampleSam",
        f"I={input_bam}",
        f"O={output_bam}",
        f"P={proportion:.6f}",
        f"RANDOM_SEED={SEED}"
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def parse_args():
    parser = argparse.ArgumentParser(description="Subsample BAM files using Picard based on target depths.")
    parser.add_argument("--bam", required=True, help="Path to BAM file.")
    parser.add_argument("--depth", nargs="+", type=float, required=True,
                        help="Target depths to simulate.")
    parser.add_argument("--out", nargs="+", required=True,
                        help="Output BAM file names for each depth.")
    parser.add_argument("--bam_size", type=int, help="Optional precomputed BAM read count.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    assert len(args.depth) == len(args.out), \
        "Number of --depth values must match number of --out outputs."

    # Determine BAM size
    BAM_SIZE = args.bam_size if args.bam_size else count_total_reads(args.bam)
    print(f"Total reads in BAM: {BAM_SIZE}")

    # Run subsampling for each depth
    for target_depth, output_bam in zip(args.depth, args.out):
        needed_reads = depth2count(target_depth)
        proportion = needed_reads / BAM_SIZE
        proportion = min(proportion, 1.0)  # cap at 1.0

        print(f"\nTarget depth: {target_depth}x")
        print(f"Estimated reads needed: {needed_reads}")
        print(f"Proportion of total: {proportion:.6f}")

        run_picard_downsample(args.bam, output_bam, proportion)

