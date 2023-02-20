from argparse import ArgumentParser
import pysam
import re

def intra_alignment_extraction(bam_file, bed_file, variant = None, extension = 50):
    # if user supplies VCF file, parse sig contents (defined below)
    if variant:
        chromosome = variant[0]
        del_start_position = variant[1]
        del_end_position = variant[1] + variant[2]

        bed = open(bed_file, "w")

        sam_file = pysam.AlignmentFile(bam_file, "rb")
        for read in sam_file.fetch(chromosome, until_eof = True):
            # get start and end position of aligned reads relative to the reference
            read_start = read.reference_start
            read_end = read_start + read.query_length

            # check for overlap between read and proposed variant in VCF file
            if max(read_start, del_start_position) <= min(read_end, del_end_position):
                current_read_position = read_start
                cigar_string = read.cigartuples

                for op, length in cigar_string:
                    if op == 2: # "D"
                        # check for overlap between CIGAR deletion and proposed variant in VCF file
                        if (max(current_read_position, del_start_position - extension) <= 
                            min(current_read_position + length, del_end_position + extension)):
                            bed.write("DEL\t" + chromosome + "\t" + str(current_read_position) + 
                                        "\t" + str(current_read_position + length) + "\n")
                        current_read_position += length

                    elif op == 0: # "M"
                        current_read_position += length

        bed.close()

def parse_vcf_file(vcf_file, chromosome):
    with open(vcf_file, "r") as file:
        contents = file.read()

        # chromosome_filter is a list of sigs of the form [(chromosome_number, start_position, length)]
        chromosome_filter = re.findall(r"\n{}.*\n".format(chromosome), contents)
        chromosome_filter = [variant.strip().split("\t") for variant in chromosome_filter]
        chromosome_filter = [(variant[0], int(variant[1]), int(re.findall(r"SVLEN=-?(\d+)", 
                                variant[-3])[0])) for variant in chromosome_filter] 
        
        return chromosome_filter

def parse_args():
    parser = ArgumentParser(description = "intra-alignment deletion signature extraction")

    parser.add_argument("-b", "--bam", default = "datasets/chr21.bam", help = "user-supplied BAM file (default: datasets/chr21.bam)")
    parser.add_argument("-c", "--chromosome", default = "chr21", help = "limits signature extraction to a particular chromosome (default: chr21)")
    parser.add_argument("-d", "--bed", default = "datasets/chr21.bed", help = "output BED file (default: chr21.bed)")
    parser.add_argument("-v", "--vcf", default = "datasets/fp.vcf", help = "user-supplied VCF file (default: datasets/fp.vcf)")

    return parser.parse_args()

def main():
    args = parse_args()

    bam_file = args.bam
    bed_file = args.bed
    chromosome = args.chromosome
    vcf_file = args.vcf

    variants = parse_vcf_file(vcf_file, chromosome)
    
    for variant in variants:
        intra_alignment_extraction(bam_file, bed_file, variant)

main()