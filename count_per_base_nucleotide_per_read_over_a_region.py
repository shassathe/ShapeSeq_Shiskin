"""
This script uses pileup function from pysam to extract reads at a given position and record the nucleotide found across
each read at the given position. The script will also record the reference nucleotide at that position. The final output
 will be a dataframe that specific the position, reference nucleotide, and the number of A,T,G,C or N nucleotides or
 insertions and deletions found across the reads mapping to that position. As input, the script accepts a bam / sam file
 , a reference fasta file ( indexed using samtools faidx ), output file location and name and the coordinates of the
 region of interest. The coordinates need to be provided in a UCSC genome browser format ( eg: chr1:1000-1040 ).
"""

# Import required packages
import pandas
import pysam
import os
import numpy
from argparse import ArgumentParser


def generate_per_base_counts(b, f, c, s, e):
    """
    The function will accept a bam file, reference fasta file and the coordinates of the region of interest. The
    function iterates over all reads found across each position in the given region and records the number of times each
     nucleotide (A,T,G,C,N) or insertions or deletions are found. The function compiles a dataframe to store the
     information and return the dataframe to be outputted.
    :param b: Full path and name of bam file as entered by the user
    :param f: Full path and name of reference fasta file as entered by the user
    :param c: Chromosome of the region of interest
    :param s: Start position of the region of interest
    :param e: End position of the region of interest
    :return: Dataframe containing pileup information across each position in the region of interest
    """
    # Read bam file and reference fasta file using pysam
    bamfile = pysam.Samfile(b)
    ref_fasta = pysam.FastaFile(f)

    # Create initial dataframe filled with zeroes
    df_temp = pandas.DataFrame(0, index=numpy.arange(s, s),
                          columns=['ref', 'A', 'T', 'G', 'C', 'del', 'ins', 'N'])

    # Use piileup to iterate over each position in the region
    for pileupcolumn in bamfile.pileup(c, s, e, stepper='all'):

        # get base position and extract reference base nucleotide
        df_temp.loc[pileupcolumn.pos, 'ref'] = ref_fasta.fetch(pileupcolumn.reference_name, int(pileupcolumn.pos) - 1,
                                                          int(pileupcolumn.pos))

        # iterate over each read at a given position and records the base nucleotide in the read. Also records if the
        # read has an insertion or deletion
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del:
                df_temp.loc[pileupcolumn.pos, 'del'] += 1
            elif pileupread.is_refskip:
                df_temp.loc[pileupcolumn.pos, 'ins'] += 1
            else:
                df_temp.loc[pileupcolumn.pos, pileupread.alignment.query_sequence[pileupread.query_position]] += 1
    return df_temp


def main():
    # Set up arguments for the script
    parser = ArgumentParser(prog='count nucleotide base along each read mapping to a given region',
                            description='Iterate over each read mapping to a given region and record the nucleotide '
                                        'base in each read at evey position within the region')
    parser.add_argument("input", help="Complete file name including the path to input bam or sam file")
    parser.add_argument("fasta", help="Complete name and path of reference fasta file. Reference fasta file must be "
                                      "indexed using samtools faidx.")
    parser.add_argument("output", help="Complete file name and path of output file")
    parser.add_argument("region", help="Region coordinates in UCSC genome browser format. eg: chr1:2000-3000")
    args = parser.parse_args()

    # Check if the input bam file exists
    if os.path.isfile(args.input):
        print('Found input bam file\n')

        # Check if the reference fasta file exists
        if os.path.isfile(args.fasta):
            print('Found reference fasta file. Proceeding with analysis...\n')

            # Extract chromosome, start position and end position of the region
            chr = args.region.split(':')[0]
            strt = int(args.region.split(':')[1].split('-')[0])
            end = int(args.region.split(':')[1].split('-')[1])

            # Analyze data through the generate_per_base_counts function
            df = generate_per_base_counts(args.input, args.fasta, chr, strt, end)
            df.to_csv(args.output, sep='\t')
        else:
            parser.error("Reference fasta file not found. Enter complete name and path of the fasta file.")
    else:
        parser.error("Input bam file not found. Enter full name and path of bam file")


if __name__ == '__main__':
    main()
