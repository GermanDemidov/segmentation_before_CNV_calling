#!/usr/bin/python3
__author__ = 'gdemidov'
import argparse
from collections import defaultdict
import os


def segment_region_into_probes(start_and_end, length_of_probe):
    list_of_coords = []
    start = int(start_and_end[0])
    end = int(start_and_end[1])

    number_of_full_probes = (end - start + 1) // length_of_probe
    if (end - start + 1) % length_of_probe > 5:
        number_of_full_probes += 1

    if number_of_full_probes == 0:
        number_of_full_probes = 1
        #print(start_and_end)
        #quit()
    if number_of_full_probes == 1:
        list_of_coords.append((start, end))
    if number_of_full_probes >= 2:
        list_of_coords.append((start, start + length_of_probe - 1))
        list_of_coords.append((end - length_of_probe, end - 1))
        num_of_probes_except_borders = number_of_full_probes - len(list_of_coords)
        if number_of_full_probes > 2:
            while float((end - length_of_probe - start)) / (num_of_probes_except_borders + 1) > length_of_probe -  10:
                # IMPORTANT STEP
                # we increase density of probes to make it looking more like real data
                num_of_probes_except_borders += 1
            if number_of_full_probes > 3:
                while float((end - length_of_probe - start)) / (num_of_probes_except_borders + 1) > round(2/3 * length_of_probe):
                    num_of_probes_except_borders += 1
                if num_of_probes_except_borders > 2 and float((end - length_of_probe - start)) / (num_of_probes_except_borders + 1) < round(2/3 * length_of_probe) - 1:
                    num_of_probes_except_borders -= 1
        if num_of_probes_except_borders > 0:
            length_of_difference = float((end - length_of_probe - start)) / (num_of_probes_except_borders + 1)
            for i in range(num_of_probes_except_borders):
                list_of_coords.append((start + round((i + 1) * length_of_difference), start + round((i + 1) * length_of_difference) + length_of_probe - 1))


    list_of_coords = sorted(list_of_coords)
    return list_of_coords


def make_probs_from_bed(bed_file, probes_file, windowLength):
    """
    We take bed files if we do not have exact probes coordintes and make it probe like (of fixed or approx fixed length)
    and extract these "emulated" probes from the reference genome
    :param bed_file: bed files - three columns, first - chrom, second - start, third - end
    :param probes_file: output file
    :param windowLength: length of probe used for hybridization
    :return: we just write everything into the file
    """
    probes = defaultdict(list)
    with open(bed_file, "r") as f:
        for line in f:
            first_symbol = line.split("\t")
            if first_symbol[0].isdigit() or first_symbol[0] in ("X","Y"):
                line = "chr" + line
            if line.startswith("chr"):
                splitted_string = line.split()
                chrom = splitted_string[0]
                coords_and_annotation = splitted_string[1:3]
                #coords_and_annotation.append(splitted_string[column_with_anotation].split(",")[0])
                probes[chrom].append(coords_and_annotation)


    dict_of_generated_probes = defaultdict(list)
    counter_of_probes = 0
    counter_of_regions = 0

    output_bed_for_coverage_calculation = defaultdict(list)
    with open(probes_file, "w") as f:
        for key_chrom in sorted(probes.keys()):
            print(key_chrom)
            with open(probes_file, "a") as f:
                for elem in probes[key_chrom]:

                    counter_of_regions += 1
                    emulated_probes = segment_region_into_probes(elem[:2], windowLength)
                    starts_and_ends = []
                    starts_and_ends_prepared = []
                    for probe in emulated_probes:
                        dict_of_generated_probes[key_chrom].append((int(probe[0]), int(probe[1])))
                        starts_and_ends.append(int(probe[0]))
                        starts_and_ends.append(int(probe[1]))
                        f.write("\t".join([key_chrom, str(probe[0]), str(probe[1])]) + "\n")
                        counter_of_probes += 1
                        if counter_of_probes % 100000 == 0:
                            print(counter_of_probes)
                    sorted_starts_ends = sorted(starts_and_ends)
                    if len(sorted_starts_ends) > 2:
                        for i in range(len(sorted_starts_ends) - 1):
                            if (sorted_starts_ends[i + 1] - sorted_starts_ends[i] > 1):
                                starts_and_ends_prepared.append([sorted_starts_ends[i], sorted_starts_ends[i + 1]])
                    else:
                        starts_and_ends_prepared.append(sorted_starts_ends)
                    #for elem in starts_and_ends_prepared:
                    #    if elem[1] - elem[0] < 2:
                    #        print(elem)
                    output_bed_for_coverage_calculation[key_chrom].extend(starts_and_ends_prepared)

            if counter_of_probes > 10**40:
                break
    probes_name = probes_file.split(".")
    for_coverage_name = ".".join(probes_name[:-1]) + ".for_coverage.bed"
    with open(for_coverage_name, "w") as f:
        for key in output_bed_for_coverage_calculation:
            for elem in output_bed_for_coverage_calculation[key]:
                f.write(key + "\t" + str(elem[0]) + "\t" + str(elem[1]) + "\n")











def main():
    parser = argparse.ArgumentParser(description='Process command line for tiling.py.')
    parser.add_argument('--bed', '-b',
                        default='regions.txt', type=str, required=True, dest='regions_file',
                        help="Bed file with regions")
    parser.add_argument('--output', '-o',
                        default='probes_from_bed.txt', type=str, required=False, dest='probes_file',
                        help="Output prboes file")
    parser.add_argument('--probLen', '-pl',
                        default="120", type=str, required=False, dest='probLen',
                        help="Length of probe used for enrichment")

    args = parser.parse_args()
    bed_file = args.regions_file
    if not args.probes_file == 'probes_from_bed.txt':
        probes_file = args.probes_file
    else:
        probes_file = ".".join(args.regions_file.split(".")[:-1]) + ".probes.bed"
    windowLength = int(args.probLen)
    print(probes_file)
    make_probs_from_bed(bed_file, probes_file, windowLength)


main()
