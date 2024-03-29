import argparse
import matplotlib
import matplotlib.image as img
import matplotlib.pyplot as plt
from multiprocessing import Process
import numpy as np
import os
import pandas as pd
import pysam
import re

matplotlib.use("agg")

def parse_vcf_file(vcf_file, chromosome):
    with open(vcf_file, "r") as file:
        contents = file.readlines()

        # chromosome_filter is a list of sigs of the form [(chromosome_number, start_position, length)]
        chromosome_filter = [line for line in contents if line[0 : line.find("\t")] == chromosome]
        chromosome_filter = [variant.strip().split("\t") for variant in chromosome_filter]
        chromosome_filter = [(variant[0], int(variant[1]), int(re.findall(r"SVLEN=-*?(\d+)", 
                              variant[-3])[0])) for variant in chromosome_filter]
        
        return chromosome_filter

def fetch_split_alignments(variant_type, chromosome, cache, duplication_factor = 2, variant_cap = 100000):
    duplicate_reads = {}
    signatures = []

    sam_file = cache[0]
    indexed_sam_file = cache[1]
    
    # fetch duplicate reads in BAM / SAM file
    for read in sam_file.fetch(chromosome, until_eof = True):
        duplicate_reads[read.query_name] = duplicate_reads.get(read.query_name, 0) + 1

    for read_name, count in duplicate_reads.items():
        if 1 < count <= duplication_factor:
            reads = indexed_sam_file.find(read_name)

            points_of_interest = []

            for read in reads:
                bias = 0
                for cigar in read.cigartuples:
                    if (cigar[0] == 0 or cigar[0] == 2 or 
                        cigar[0] == 7 or cigar[0] == 8): # M or D or = or X
                        bias += cigar[1]

                orientation = "+" if read.is_forward else "-"

                points_of_interest.append((read.reference_start, read.reference_start + bias, read.query_alignment_start, 
                                           read.query_alignment_end, orientation, read.reference_name))
            
            for index in range(len(points_of_interest) - 1):
                first_segment = points_of_interest[index]
                second_segment = points_of_interest[index + 1]
                
                # do the reads have the same orientation / are on the same chromosome?
                if first_segment[4] == second_segment[4] and first_segment[5] == second_segment[5] == chromosome:
                    # utilize CuteSV heuristic for identifying deletion / insertion signatures
                    difference_distance = (second_segment[0] - first_segment[1]) - (second_segment[2] - first_segment[3])
                    difference_overlap = first_segment[1] - second_segment[0]

                    # compose deletion signature as identified by split alignment
                    if variant_type == "DEL":
                        if difference_overlap < 30 and 30 <= difference_distance <= variant_cap:
                            signatures.append((first_segment[1], difference_distance, read_name))
                    
                    # compose insertion signature as identified by split alignment
                    elif variant_type == "INS":
                        if difference_overlap < 30 and -1 * variant_cap <= difference_distance <= -30:
                            signatures.append(((first_segment[1] + second_segment[0]) // 2, -1 * difference_distance, read_name))

                    else:
                        raise ValueError("unsupported variant type provided")
    
    return signatures

def intra_alignment_extraction(variant_type, bam_file, bed, variant = None, extension = 50):
    # if user supplies VCF file, parse sig contents (defined below)
    if variant:
        chromosome = variant[0]
        variant_start_position = variant[1]
        variant_end_position = variant[1] + variant[2]

        bed_file = open(bed + "/{}_{}_{}.bed".format(chromosome, variant_start_position, variant_end_position), "w")
        bed_file.write("{}\t{}\t{}\t{}\t{}\n".format("CHROMOSOME", "START", "END", "READ", "TYPE"))

        sam_file = pysam.AlignmentFile(bam_file, "rb")

        for read in sam_file.fetch(chromosome, until_eof = True):
            # get start and end position of aligned reads relative to the reference
            read_start = read.reference_start
            read_end = read_start + read.query_length

            # check for overlap between read and proposed variant in VCF file
            if max(read_start, variant_start_position) <= min(read_end, variant_end_position):
                current_read_position = read_start
                cigar_string = read.cigartuples

                for op, length in cigar_string:
                    if variant_type == "DEL":
                        if op == 0: # "M"
                            current_read_position += length

                        elif op == 2: # "D"
                            # check for overlap between CIGAR deletion and proposed variant in VCF file
                            if (max(current_read_position, variant_start_position - extension) <= 
                                min(current_read_position + length, variant_end_position + extension)):
                                bed_file.write(chromosome + "\t" + str(current_read_position) + "\t" + str(current_read_position 
                                               + length) + "\t" + read.query_name + "\t" + "INTRA_DEL" + "\n")
                                
                            current_read_position += length                            
                    
                    elif variant_type == "INS":
                        if op == 0: # "M"
                            current_read_position += length
                        
                        elif op == 1: # "I"
                            if (max(current_read_position, variant_start_position - extension) <= 
                                min(current_read_position + length, variant_end_position + extension)):
                                bed_file.write(chromosome + "\t" + str(current_read_position) + "\t" + str(current_read_position 
                                               + length) + "\t" + read.query_name + "\t" + "INTRA_INS" + "\n")
                            
                            current_read_position += length

                        elif op == 2: # "D"
                            current_read_position += length

                    else:
                        raise ValueError("unsupported variant type provided")

        bed_file.close()
    
    else:
        # TODO: add logic to handle this case by parsing BAM directly
        raise FileNotFoundError("VCF file was not supplied")

def inter_alignment_extraction(variant_type, split_read_signatures, bed, variant = None, extension = 50):
    chromosome = variant[0]
    variant_start_position = variant[1]
    variant_end_position = variant[1] + variant[2]

    bed_file = open(bed + "/{}_{}_{}.bed".format(chromosome, variant_start_position, variant_end_position), "a")

    for signature in split_read_signatures:
        signature_interval = (signature[0], signature[0] + signature[1])

        # determine whether variant and signature intervals overlap
        if (max(variant_start_position - extension, signature_interval[0]) <= 
            min(variant_end_position + extension, signature_interval[1])):
            bed_file.write(chromosome + "\t" + str(signature_interval[0]) + "\t" + str(signature_interval[1]) + 
                           "\t" + signature[2] + "\t" + "INTER_{}".format(variant_type) + "\n")
    
    bed_file.close()

def visualize_alignments(variant_type, images, bed, variant):
    chromosome = variant[0]
    variant_start_position = variant[1]
    variant_end_position = variant[1] + variant[2]
    
    read_alignments = pd.read_csv(bed + "/{}_{}_{}.bed".format(chromosome, variant_start_position, 
                                                               variant_end_position), sep = "\t")
    read_alignments["HEIGHT"] = read_alignments.groupby("READ").ngroup() + 2

    plt.plot([variant_start_position, variant_end_position], [1, 1], color = "blue")
    for _, row in read_alignments.iterrows():
        plt.plot([row["START"], row["END"]], [row["HEIGHT"], row["HEIGHT"]], 
                 color = "orange" if row["TYPE"] == "INTRA_{}".format(variant_type) else "green")
    
    plt.xticks(rotation = "vertical")
    plt.tick_params(left = False, labelleft = False)
    plt.ticklabel_format(axis = "x", style = "sci", useOffset = False)
    
    plt.tight_layout()
    
    plt.savefig(images + "/signatures/{}_{}_{}.png".format(chromosome, variant_start_position, 
                                                           variant_end_position))
    plt.clf()

def encode_variant_as_matrix(variant_type, bed, variant, normalize_by_padding = False, 
                             normalized_width = 10000, normalized_height = 60, extension = 50):
    variant_matrix = {}

    chromosome = variant[0]
    variant_start_position = variant[1]
    variant_end_position = variant[1] + variant[2]

    def vectorize(read_start, read_end, extraction):
        read_vector = [0] * ((variant_end_position - variant_start_position) + 2 * extension)

        start_mask = max(variant_start_position - extension, read_start) - (variant_start_position - extension)
        end_mask = min(variant_end_position + extension, read_end) - (variant_start_position - extension)

        mask = ([1] * (end_mask - start_mask) if extraction == "INTRA_{}".
                format(variant_type) else [2] * (end_mask - start_mask))
        
        read_vector[start_mask : end_mask] = mask
        return read_vector

    with open(bed + "/{}_{}_{}.bed".format(chromosome, variant_start_position, 
                                           variant_end_position), "r") as bed_file:
        _ = bed_file.readline() # ignore BED header information

        for read in bed_file:
            read = read.split("\t")
            read_vector = vectorize(int(read[1]), int(read[2]), read[4].rstrip())
            np.bitwise_or(variant_matrix.setdefault(read[3], np.array(read_vector)), 
                          np.array(read_vector), out = variant_matrix[read[3]])
    
    matrix = np.matrix(list(variant_matrix.values()))

    if normalize_by_padding:
        unnormalized_rows = np.shape(matrix)[0]
        unnormalized_columns = np.shape(matrix)[1]

        total_vertical_padding = normalized_height - unnormalized_rows
        vertical_padding_per_side = total_vertical_padding // 2
        vertical_padding_per_side = ((vertical_padding_per_side, vertical_padding_per_side + 1) if total_vertical_padding % 2 != 0 else
                                     (vertical_padding_per_side, vertical_padding_per_side))

        total_horizontal_padding = normalized_width - unnormalized_columns
        horizontal_padding_per_side = total_horizontal_padding // 2
        horizontal_padding_per_side = ((horizontal_padding_per_side, horizontal_padding_per_side + 1) if total_horizontal_padding % 2 != 0 else
                                       (horizontal_padding_per_side, horizontal_padding_per_side))
        
        if total_horizontal_padding >= 0 and total_vertical_padding >= 0:    
            matrix = np.pad(matrix, (vertical_padding_per_side, horizontal_padding_per_side))

        elif total_horizontal_padding >= 0:
            matrix = np.pad(matrix, ((0, 0), horizontal_padding_per_side))
            matrix = matrix[abs(vertical_padding_per_side[1]) : 
                            abs(vertical_padding_per_side[1]) + normalized_height, :]

        elif total_vertical_padding >= 0:
            matrix = np.pad(matrix, (vertical_padding_per_side, (0, 0)))
            matrix = matrix[:, abs(horizontal_padding_per_side[1]) : 
                            abs(horizontal_padding_per_side[1]) + normalized_width]

        else:
            matrix = matrix[abs(vertical_padding_per_side[1]) : 
                            abs(vertical_padding_per_side[1]) + normalized_height, :]
            matrix = matrix[:, abs(horizontal_padding_per_side[1]) : 
                            abs(horizontal_padding_per_side[1]) + normalized_width]

        assert(np.shape(matrix)[0] == normalized_height)
        assert(np.shape(matrix)[1] == normalized_width)

    return matrix

def generate_encoding(images, variant_matrix, variant, encoding_format):
    chromosome = variant[0]
    variant_start_position = variant[1]
    variant_end_position = variant[1] + variant[2]

    if encoding_format == "plot":
        image = np.zeros(tuple(list(variant_matrix.shape) + [3]), dtype = np.uint8)

        for read_number, read in enumerate(variant_matrix):
            for bp_number, bp in enumerate(read.flat):
                if bp != 0:
                    image[read_number][bp_number] = [int((bp != 2) * 255), int((bp != 1) * 255), 0]

        img.imsave(images + "/matrices/{}_{}_{}.png".format(chromosome, variant_start_position, 
                                                            variant_end_position), image)
    
    else:
        intra_encoding_counts = np.sum(np.where(variant_matrix == 1, variant_matrix, 0), axis = 0)
        inter_encoding_counts = np.sum(np.where(variant_matrix == 2, variant_matrix, 0), axis = 0)

        image = np.zeros((max(max(intra_encoding_counts), max(inter_encoding_counts)) + 1, 
                          variant_matrix.shape[1], 3), dtype = np.uint8)
        
        intra_encoding_counts = list(enumerate(intra_encoding_counts))
        inter_encoding_counts = list(enumerate(inter_encoding_counts))

        for column in range(variant_matrix.shape[1]):
            intra_coordinate = intra_encoding_counts[column]
            inter_coordinate = inter_encoding_counts[column]

            if intra_coordinate[1] != 0:
                image[intra_coordinate[1]][intra_coordinate[0]][0] = 255
            
            if inter_coordinate[1] != 0:
                image[inter_coordinate[1]][inter_coordinate[0]][1] = 255
        
        image = np.flipud(image)

        img.imsave(images + "/matrices/{}_{}_{}.png".format(chromosome, variant_start_position, 
                                                            variant_end_position), image)

def parse_args():
    parser = argparse.ArgumentParser(description = "insertion & deletion signature extraction")

    parser.add_argument("-b", "--bam", default = "data/chr21.bam", help = "user-supplied BAM file; use the keyword 'all' for the entire genome (default: data/chr21.bam)")
    parser.add_argument("-c", "--chromosomes", default = "chr21", help = "limits signature extraction to particular chromosomes; \
                                                                          specify as a comma separated list or using the keyword 'all' for the entire genome (default: chr21)")
    parser.add_argument("-d", "--bed", default = "data/bed", help = "output BED directory (default: data/bed)")
    parser.add_argument("-f", "--format", default = "counts", choices = ["counts, plot"], help = "model input options (default: counts)")
    parser.add_argument("-i", "--images", default = "data/images", help = "output image directory (default: data/images)")
    parser.add_argument("-n", "--normalize", action = "store_true", help = "normalize matrices to a fixed width and height")
    parser.add_argument("-t", "--type", default = "DEL", choices = ["DEL", "INS"], help = "structural variant type (default: DEL)")
    parser.add_argument("-v", "--vcf", default = "data/fp.vcf", help = "user-supplied VCF file (default: data/fp.vcf)")

    return parser.parse_args()

def launch_chromosome_extraction(bam_file, chromosome, base_bed, base_images, 
                                 variant_type, vcf_file, normalize, encoding_format = "counts"):
    sam_file = pysam.AlignmentFile(bam_file, "rb")
    indexed_sam_file = pysam.IndexedReads(sam_file)
    indexed_sam_file.build()

    cache = [sam_file, indexed_sam_file]

    bed = os.path.join(base_bed, chromosome)
    images = os.path.join(base_images, chromosome)

    os.makedirs(bed, exist_ok = True)
    os.makedirs(os.path.join(images, "matrices"), exist_ok = True)
    os.makedirs(os.path.join(images, "signatures"), exist_ok = True)

    variants = parse_vcf_file(vcf_file, chromosome)
    split_read_signatures = fetch_split_alignments(variant_type, chromosome, cache)
    
    for variant in variants:
        print("analyzing variant on {} with start position {} and end position {}".format(variant[0], variant[1], 
                                                                                          variant[1] + variant[2]))
        
        intra_alignment_extraction(variant_type, bam_file, bed, variant)
        inter_alignment_extraction(variant_type, split_read_signatures, bed, variant)
        visualize_alignments(variant_type, images, bed, variant)
        generate_encoding(images, encode_variant_as_matrix(variant_type, bed, variant, normalize), variant, encoding_format)

def main():
    args = parse_args()

    bam_file = args.bam
    chromosomes = (["chr{}".format(chromosome_number) for chromosome_number in range(1, 23)] 
                        if args.chromosomes == "all" else args.chromosomes.split(","))
    base_bed = args.bed
    encoding_format = args.format
    base_images = args.images
    variant_type = args.type
    vcf_file = args.vcf
    normalize = args.normalize

    chromosome_processes = []
    for chromosome in chromosomes:
        chromosome_processes.append(Process(target = launch_chromosome_extraction, args = (bam_file, chromosome, base_bed, base_images, 
                                                                                           variant_type, vcf_file, normalize, encoding_format)))
        chromosome_processes[-1].start()
    
    for chromosome_process in chromosome_processes:
        chromosome_process.join()

if __name__ == "__main__":
    main()