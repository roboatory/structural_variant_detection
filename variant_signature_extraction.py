from argparse import ArgumentParser
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import re

def parse_vcf_file(vcf_file, chromosome):
    with open(vcf_file, "r") as file:
        contents = file.readlines()

        # chromosome_filter is a list of sigs of the form [(chromosome_number, start_position, length)]
        chromosome_filter = [line for line in contents if line.startswith(chromosome)]
        chromosome_filter = [variant.strip().split("\t") for variant in chromosome_filter]
        chromosome_filter = [(variant[0], int(variant[1]), int(re.findall(r"SVLEN=-?(\d+)", 
                              variant[-3])[0])) for variant in chromosome_filter]
        
        return chromosome_filter

def fetch_split_alignments(bam_file, chromosome, duplication_factor = 2, variant_cap = 100000):
    duplicate_reads = {}
    deletion_signatures = []

    sam_file = pysam.AlignmentFile(bam_file, "rb")
    indexed_sam_file = pysam.IndexedReads(sam_file)
    indexed_sam_file.build()

    # fetch duplicate reads in BAM / SAM file
    for read in sam_file.fetch(chromosome, until_eof = True):
        duplicate_reads[read.query_name] = duplicate_reads.get(read.query_name, 0) + 1

    for read_name, count in duplicate_reads.items():
        if 1 < count <= duplication_factor:
            reads = indexed_sam_file.find(read_name)

            points_of_interest = [(read.reference_start, read.reference_start + read.query_length, 
                                   read.query_alignment_start, read.query_alignment_end) for read in reads]
            
            for index in range(len(points_of_interest) - 1):
                first_segment = points_of_interest[index]
                second_segment = points_of_interest[index + 1]

                # utilize CuteSV heuristic for identifying deletion signatures
                difference_distance = (second_segment[0] - first_segment[1]) - (second_segment[2] - first_segment[3])
                difference_overlap = first_segment[1] - second_segment[0]

                # compose deletion signature as identified by split alignment
                if difference_overlap < 30 and 30 <= difference_distance <= variant_cap:
                    deletion_signatures.append((first_segment[1], difference_distance, read_name))
    
    return deletion_signatures

def intra_alignment_extraction(bam_file, bed, variant = None, extension = 50):
    # if user supplies VCF file, parse sig contents (defined below)
    if variant:
        chromosome = variant[0]
        del_start_position = variant[1]
        del_end_position = variant[1] + variant[2]

        bed_file = open(bed + "/{}_{}_{}.bed".format(chromosome, del_start_position, del_end_position), "w")
        bed_file.write("{}\t{}\t{}\t{}\t{}\n".format("CHROMOSOME", "START", "END", "READ", "TYPE"))

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
                            bed_file.write(chromosome + "\t" + str(current_read_position) + "\t" + str(current_read_position 
                                           + length) + "\t" + read.query_name + "\t" + "INTRA_DEL" + "\n")
                            
                        current_read_position += length

                    elif op == 0: # "M"
                        current_read_position += length

        bed_file.close()
    
    else:
        # TODO: add logic to handle this case by parsing BAM directly
        raise FileNotFoundError("VCF file was not supplied")

def inter_alignment_extraction(split_read_deletion_signatures, bed, variant = None, extension = 50):
    chromosome = variant[0]
    del_start_position = variant[1]
    del_end_position = variant[1] + variant[2]

    bed_file = open(bed + "/{}_{}_{}.bed".format(chromosome, del_start_position, del_end_position), "a")

    for deletion_signature in split_read_deletion_signatures:
        signature_interval = (deletion_signature[0], deletion_signature[0] + deletion_signature[1])

        # determine whether variant and signature intervals overlap
        if (max(del_start_position - extension, signature_interval[0]) <= 
            min(del_end_position + extension, signature_interval[1])):
            bed_file.write(chromosome + "\t" + str(signature_interval[0]) + "\t" + str(signature_interval[1]) + 
                           "\t" + deletion_signature[2] + "\t" + "INTER_DEL" + "\n")
    
    bed_file.close()

def visualize_alignments(images, bed, variant):
    chromosome = variant[0]
    del_start_position = variant[1]
    del_end_position = variant[1] + variant[2]
    
    read_alignments = pd.read_csv(bed + "/{}_{}_{}.bed".format(chromosome, del_start_position, 
                                                               del_end_position), sep = "\t")
    read_alignments["HEIGHT"] = read_alignments.groupby("READ").ngroup() + 2

    plt.plot([del_start_position, del_end_position], [1, 1], color = "blue")
    for _, row in read_alignments.iterrows():
        plt.plot([row["START"], row["END"]], [row["HEIGHT"], row["HEIGHT"]], 
                 color = "orange" if row["TYPE"] == "INTRA_DEL" else "green")
    
    plt.xticks(rotation = "vertical")
    plt.tick_params(left = False, labelleft = False)
    plt.ticklabel_format(axis = "x", style = "sci", useOffset = False)
    
    plt.tight_layout()
    
    plt.savefig(images + "/signatures/{}_{}_{}.png".format(chromosome, del_start_position, 
                                                           del_end_position))
    plt.clf()

def encode_variant_as_matrix(bed, variant, extension = 50):
    variant_matrix = {}

    chromosome = variant[0]
    del_start_position = variant[1]
    del_end_position = variant[1] + variant[2]

    def vectorize(read_start, read_end, extraction):
        read_vector = [0] * ((del_end_position - del_start_position) + 2 * extension)

        start_mask = max(del_start_position - extension, read_start) - (del_start_position - extension)
        end_mask = min(del_end_position + extension, read_end) - (del_start_position - extension)

        mask = ([1] * (end_mask - start_mask) if extraction == "INTRA_DEL" else 
                [2] * (end_mask - start_mask))
        
        read_vector[start_mask : end_mask] = mask
        return read_vector

    with open(bed + "/{}_{}_{}.bed".format(chromosome, del_start_position, 
                                           del_end_position), "r") as bed_file:
        _ = bed_file.readline() # ignore BED header information

        for read in bed_file:
            read = read.split("\t")
            read_vector = vectorize(int(read[1]), int(read[2]), read[4].rstrip())
            np.bitwise_or(variant_matrix.setdefault(read[3], np.array(read_vector)), 
                          np.array(read_vector), out = variant_matrix[read[3]])
    
    return np.matrix(list(variant_matrix.values())).astype(float)

def visualize_matrix_encoding(images, variant_matrix, variant):
    chromosome = variant[0]
    del_start_position = variant[1]
    del_end_position = variant[1] + variant[2]

    # white - no alignment, orange - intra-alignment, green - inter-alignment
    color_map = colors.ListedColormap(["white", "orange", "green"])

    plt.matshow(variant_matrix, cmap = color_map, vmin = 0., vmax = 2.)
    plt.tick_params(left = False, bottom = False, top = False, 
                    labelleft = False, labeltop = False)

    plt.savefig(images + "/matrices/{}_{}_{}.png".format(chromosome, del_start_position, 
                                                         del_end_position))
    plt.close()

def parse_args():
    parser = ArgumentParser(description = "intra-alignment deletion signature extraction")

    parser.add_argument("-b", "--bam", default = "data/chr21.bam", help = "user-supplied BAM file (default: data/chr21.bam)")
    parser.add_argument("-c", "--chromosome", default = "chr21", help = "limits signature extraction to a particular chromosome (default: chr21)")
    parser.add_argument("-d", "--bed", default = "data/bed", help = "output BED directory (default: data/bed)")
    parser.add_argument("-i", "--images", default = "data/images", help = "output image directory (default: data/images")
    parser.add_argument("-v", "--vcf", default = "data/fp.vcf", help = "user-supplied VCF file (default: data/fp.vcf)")

    return parser.parse_args()

def main():
    args = parse_args()

    bam_file = args.bam
    bed = args.bed
    chromosome = args.chromosome
    vcf_file = args.vcf
    images = args.images

    variants = parse_vcf_file(vcf_file, chromosome)
    split_read_deletion_signatures = fetch_split_alignments(bam_file, chromosome)
    
    for variant in variants:
        intra_alignment_extraction(bam_file, bed, variant)
        inter_alignment_extraction(split_read_deletion_signatures, bed, variant)
        visualize_alignments(images, bed, variant)
        visualize_matrix_encoding(images, encode_variant_as_matrix(bed, variant), variant)

main()