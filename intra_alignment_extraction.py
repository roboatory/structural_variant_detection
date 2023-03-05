from argparse import ArgumentParser
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import re

def visualize_alignments(bed, variant, images):
    chromosome = variant[0]
    del_start_position = variant[1]
    del_end_position = variant[1] + variant[2]
    
    read_alignments = pd.read_csv(bed + "/{}_{}_{}.bed".format(chromosome, 
                          del_start_position, del_end_position), sep = "\t")
    read_alignments["HEIGHT"] = read_alignments.groupby("READ").ngroup() + 2

    plt.plot([del_start_position, del_end_position], [1, 1], color = "blue")
    for _, row in read_alignments.iterrows():
        plt.plot([row["START"], row["END"]], [row["HEIGHT"], row["HEIGHT"]], color = "orange")

    # TODO: make figure nicer
    
    plt.savefig(images + "/{}_{}_{}.png".format(chromosome, del_start_position, del_end_position))

def visualize_matrix(matrix, variant, images):
    chromosome = variant[0]
    del_start_position = variant[1]
    del_end_position = variant[1] + variant[2]

    # TODO: make extensible colormap

    color_map = ListedColormap(["r", "g"])
    plt.matshow(matrix, cmap = color_map)
    plt.savefig(images + "/{}_{}_{}_matrix.png".format(chromosome, del_start_position, del_end_position))

def vectorize(variant, extension_length, read_sv_start, read_sv_length, read_vector):
    del_start_position = variant[1]
    del_end_position = variant[1] + variant[2]

    read_sv_end = read_sv_start + read_sv_length

    read_sv_rep_start = max(del_start_position - extension_length, read_sv_start) - (del_start_position - extension_length)
    read_sv_rep_end = min(del_end_position + extension_length, read_sv_end) - (del_start_position - extension_length)

    mask = [1] * (read_sv_rep_end - read_sv_rep_start)
    read_vector[read_sv_rep_start : read_sv_rep_end] = mask

    return True

def intra_alignment_extraction(bam_file, bed, variant = None, extension = 50):
    # if user supplies VCF file, parse sig contents (defined below)
    if variant:
        chromosome = variant[0]
        del_start_position = variant[1]
        del_end_position = variant[1] + variant[2]

        bed_file = open(bed + "/{}_{}_{}.bed".format(chromosome, del_start_position, del_end_position), "w")
        bed_file.write("{}\t{}\t{}\t{}\t{}\n".format("CHROMOSOME", "START", "END", "READ", "TYPE"))

        sam_file = pysam.AlignmentFile(bam_file, "rb")
        image_matrices = []

        for read in sam_file.fetch(chromosome, until_eof = True):
            # get start and end position of aligned reads relative to the reference
            read_start = read.reference_start
            read_end = read_start + read.query_length

            # check for overlap between read and proposed variant in VCF file
            if max(read_start, del_start_position) <= min(read_end, del_end_position):
                current_read_position = read_start
                cigar_string = read.cigartuples

                read_vector = [0] * ((del_end_position - del_start_position) + 2 * extension)
                encountered_overlap = False

                for op, length in cigar_string:
                    if op == 2: # "D"
                        # check for overlap between CIGAR deletion and proposed variant in VCF file
                        if (max(current_read_position, del_start_position - extension) <= 
                            min(current_read_position + length, del_end_position + extension)):
                            encountered_overlap = vectorize(variant, extension, current_read_position, length, read_vector)
                            bed_file.write(chromosome + "\t" + str(current_read_position) + "\t" + str(current_read_position 
                                               + length) + "\t" + read.query_name + "\t" + "DEL" + "\n")
                        current_read_position += length

                    elif op == 0: # "M"
                        current_read_position += length

                if encountered_overlap:       
                    image_matrices.append(read_vector)

        bed_file.close()
        return image_matrices

def inter_alignment_extraction(bam_file, chromosome, duplication_factor = 2):
    duplicate_reads = {}

    sam_file = pysam.AlignmentFile(bam_file, "rb")
    indexed_sam_file = pysam.IndexedReads(sam_file)
    indexed_sam_file.build()

    # fetch duplicate reads in BAM / SAM file
    for read in sam_file.fetch(chromosome, until_eof = True):
        duplicate_reads[read.query_name] = duplicate_reads.get(read.query_name, 0) + 1

    for read_name, count in duplicate_reads.items():
        if 1 < count <= duplication_factor:
            reads = indexed_sam_file.find(read_name)

            # extract metadata about each read's alignment
            points_of_interest = [(read.reference_start, read.reference_start + read.query_length, 
                                   read.query_alignment_start, read.query_alignment_end) for read in reads]
            print(points_of_interest)

            for index in range(len(points_of_interest) - 1):
                first_segment = points_of_interest[index]
                second_segment = points_of_interest[index + 1]

                # utilize CuteSV heuristic for identifying deletion signatures
                difference_distance = (second_segment[0] - first_segment[1]) - (second_segment[2] - first_segment[3])
                difference_overlap = first_segment[1] - second_segment[0]

                if difference_overlap < 30 and difference_distance >= 30:
                    print(first_segment[1], difference_distance, read_name)

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
    parser.add_argument("-d", "--bed", default = "datasets/bed", help = "output BED directory (default: datasets/bed)")
    parser.add_argument("-i", "--images", default = "datasets/images", help = "output image directory (default: datasetes/images")
    parser.add_argument("-v", "--vcf", default = "datasets/fp.vcf", help = "user-supplied VCF file (default: datasets/fp.vcf)")

    return parser.parse_args()

def main():
    args = parse_args()

    bam_file = args.bam
    bed = args.bed
    chromosome = args.chromosome
    vcf_file = args.vcf
    images = args.images

    # variants = parse_vcf_file(vcf_file, chromosome)
    
    # for variant in variants:
    #     intra_variant_matrix = intra_alignment_extraction(bam_file, bed, variant)
    #     visualize_alignments(bed, variant, images)
    #     visualize_matrix(intra_variant_matrix, variant, images)

    inter_alignment_extraction(bam_file, chromosome)

main()