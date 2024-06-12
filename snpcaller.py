import pysam
import argparse
import pandas as pd
from math import log, exp


def compute_posterior(ref_probs, alt_probs, allele_freq):
    """Compute posterior probabilities."""
    prob_AA = (1 - allele_freq)**2
    prob_BB = allele_freq**2
    prob_AB = 2 * allele_freq * (1 - allele_freq)

    likelihood_AA = sum(log(1 - e) for e in ref_probs) + sum(log(e) for e in alt_probs)
    likelihood_BB = sum(log(e) for e in ref_probs) + sum(log(1 - e) for e in alt_probs)
    likelihood_AB = len(ref_probs + alt_probs) * log(0.5)

    denominator = exp(likelihood_AA + log(prob_AA)) + exp(likelihood_BB + log(prob_BB)) + exp(likelihood_AB + log(prob_AB))

    post_AA = exp(likelihood_AA + log(prob_AA)) / denominator
    post_AB = exp(likelihood_AB + log(prob_AB)) / denominator
    post_BB = exp(likelihood_BB + log(prob_BB)) / denominator

    return post_AA, post_AB, post_BB


def gather_reads(bam, chrom, pos, ref_base, alt_base):
    """Retrieve reads at SNP position."""
    nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    ref_error_probs = []
    alt_error_probs = []

    for pileup_column in bam.pileup(chrom, pos-1, pos, min_base_quality=0):
        if pileup_column.pos == pos-1:
            for pileup_read in pileup_column.pileups:
                if not pileup_read.is_del and not pileup_read.is_refskip:
                    base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                    quality = pileup_read.alignment.query_qualities[pileup_read.query_position]
                    error_prob = 10 ** (-quality / 10.0)

                    if base in nucleotide_counts:
                        nucleotide_counts[base] += 1
                        if base == ref_base.upper():
                            ref_error_probs.append(error_prob)
                        elif base == alt_base.upper():
                            alt_error_probs.append(error_prob)

    return ref_error_probs, alt_error_probs, sum(nucleotide_counts.values())


def filter_snps(snp_file, meta_file):
    """Load and filter SNPs using metadata."""
    snps_data = pd.read_csv(snp_file, sep='\t')
    meta_data = pd.read_csv(meta_file, sep='\t', names=["chr", "pos", "ref"], skiprows=1)

    meta_data['pos'] = meta_data['pos'].astype(int)
    snps_data['pos'] = snps_data['pos'].astype(int)

    filtered_snps = snps_data[snps_data['pos'].isin(meta_data['pos'])]

    return filtered_snps



def analyze_snps(bam_path, snps_data):
    """Process each SNP and compute probabilities."""
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    analysis_results = []

    for _, snp in snps_data.iterrows():
        ref_errors, alt_errors, total_bases = gather_reads(bam_file, str(snp['chr']), snp['pos'], snp['ref'], snp['alt'])
        post_ref, post_het, post_alt = compute_posterior(ref_errors, alt_errors, snp['maf'])

        if post_ref > post_alt:
            predicted_genotype = snp['ref'] * 2
            max_posterior = post_ref
        else:
            predicted_genotype = snp['alt'] * 2
            max_posterior = post_alt

        analysis_results.append([snp['chr'], snp['pos'], predicted_genotype, max_posterior, total_bases])

    return analysis_results



def save_to_csv(results):
    """Save results to CSV."""
    output_file = "output_results.csv"
    results_df = pd.DataFrame(results, columns=['chromosome', 'position', 'putative_genotype', 'posterior_probability', 'n_reads'])
    results_df.to_csv(output_file, index=False)
    print("Results saved to CSV.")


def index_bam(bam_file):
    try:
        pysam.index(bam_file)
        print("Indexing completed successfully.")
    except Exception as e:
        print(f"Error indexing BAM file: {e}")
        exit(1)


def run_snp_caller(bam_path, metadata_path):
    """Run the SNP caller pipeline."""
    snp_file = 'putatative_snps.tsv'

    index_bam(bam_path)
    filtered_snps_data = filter_snps(snp_file, metadata_path)
    results = analyze_snps(bam_path, filtered_snps_data)
    save_to_csv(results)

def main():
    parser = argparse.ArgumentParser(description="SNP calling from BAM file using metadata.")
    parser.add_argument("bam_file", help="Input BAM file.")
    parser.add_argument("metadata_file", help="Metadata file for SNP filtering.")

    args = parser.parse_args()

    run_snp_caller(args.bam_file, args.metadata_file)

if __name__ == "__main__":
    main()