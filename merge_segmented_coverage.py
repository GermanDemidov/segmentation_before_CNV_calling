__author__ = 'gdemidov'
import argparse
from collections import defaultdict
import operator


def intersection(interval1, interval2):
    #s1 = interval1[0]
    #s2 = interval2[0]
    #e1 = interval1[1]
    #e2 = interval2[1]
    m = min(interval1[1], interval2[1])
    s = max(interval1[0], interval2[0])
    #intersected_part = [max(s1, s2), min(e1, e2)]
    if m > s:
        return [s, m]
    else:
        return False

def divide_coverage(probes_from_bed, coverages_from_segmented):
    #print()
    #print(probes_from_bed)
    #print(coverages_from_segmented)
    if len(probes_from_bed) == len(coverages_from_segmented):
        return([coverages_from_segmented[i][2] for i in range(len(probes_from_bed))])
    coverages = [0 for i in range(len(probes_from_bed))]
    for i in range(len(probes_from_bed) - 1):
        #print(probes_from_bed[i])
        #print(probes_from_bed[i + 1])
        first_cov = coverages_from_segmented[2 * i][2] * (coverages_from_segmented[2 * i][1] - coverages_from_segmented[2 * i][0])
        second_cov = coverages_from_segmented[2 * i + 2][2] * (coverages_from_segmented[2 * i + 2][1] - coverages_from_segmented[2 * i + 2][0])
        middle_cov = coverages_from_segmented[2 * i + 1][2] * (coverages_from_segmented[2 * i + 1][1] - coverages_from_segmented[2 * i + 1][0])
        #print(first_cov)
        #print(second_cov)
        #print(middle_cov)
        if (first_cov > 0):
            if coverages[i] == 0:
                coverages[i] += first_cov
            coverages[i] += first_cov * middle_cov / (first_cov + second_cov)
        if (second_cov > 0):
            if coverages[i + 1] == 0:
                coverages[i + 1] += second_cov
            coverages[i + 1] += second_cov * middle_cov / (first_cov + second_cov)
    for i in range(len(probes_from_bed)):
        coverages[i] /= (probes_from_bed[i][1] - probes_from_bed[i][0])
    #print(coverages_from_segmented)
    #print(coverages)
    return(coverages)



def merge_coverage_file(bed_file, coverage_file, output_file):
    coverages = defaultdict(list)
    bedFile = defaultdict(list)
    header = ""
    with open(coverage_file) as f:
        header = f.readline().strip()
        for line in f:
            splitted_cov = line.strip().split("\t")
            splitted_cov[1] = int(splitted_cov[1])
            splitted_cov[2] = int(splitted_cov[2])
            splitted_cov[3] = float(splitted_cov[3])
            coverages[splitted_cov[0]].append(splitted_cov[1:])
    with open(bed_file) as f:
         for line in f:
            splitted_cov = line.strip().split("\t")
            splitted_cov[1] = int(splitted_cov[1])
            splitted_cov[2] = int(splitted_cov[2])
            bedFile[splitted_cov[0]].append(splitted_cov[1:])
    for key in coverages.keys():
        coverages[key] = sorted(coverages[key], key=operator.itemgetter(0, 1))
        bedFile[key] = sorted(bedFile[key], key=operator.itemgetter(0, 1))
    lines_to_write = [header]
    for chrom in coverages.keys():
        coverages_location = 0
        cluster = [bedFile[chrom][0]]
        bedFileIncluded = [False for i in range(len(bedFile[chrom]))]
        bedFileIncluded[0] = True
        for i in range(1, len(bedFile[chrom])):
            if intersection(bedFile[chrom][i], cluster[-1]):
                cluster.append(bedFile[chrom][i])
                bedFileIncluded[i] = True
            if (not intersection(bedFile[chrom][i], cluster[-1])) or (i == len(bedFile[chrom]) - 1):
                min_elem = cluster[0][0]
                max_elem = cluster[-1][1]
                coveragesFromCluster = []
                coverages_location_tmp = coverages_location
                for j in range(coverages_location, len(coverages[chrom])):
                    if intersection((min_elem, max_elem), coverages[chrom][j]):
                        coveragesFromCluster.append(coverages[chrom][j])
                        coverages_location_tmp = j
                    if max_elem < coverages[chrom][j][0]:
                        break
                coverages_location = coverages_location_tmp
                    #else:
                    #    coverages_location -= 1
                    #    coverages_location = max(0, coverages_location)
                    #    break
                divided_coverage = divide_coverage(cluster, coveragesFromCluster)
                for l in range(len(divided_coverage)):
                    lines_to_write.append("\t".join([chrom, str(cluster[l][0]), str(cluster[l][1]), str(round(divided_coverage[l], 3))]))
                if i == len(bedFile[chrom]) - 1:
                    if cluster[-1][0] != bedFile[chrom][i][0] or cluster[-1][1] != bedFile[chrom][i][1]:
                        cluster = [bedFile[chrom][i]]
                        coveragesFromCluster = [coverages[chrom][-1]]
                        divided_coverage = divide_coverage(cluster, coveragesFromCluster)
                        lines_to_write.append("\t".join([chrom, str(cluster[0][0]), str(cluster[0][1]), str(round(divided_coverage[0], 3))]))
                cluster = [bedFile[chrom][i]]
                bedFileIncluded[i] = True

    with open(output_file, "w") as f:
        for line in lines_to_write:
            f.write(line + "\n")





def main():
    parser = argparse.ArgumentParser(description='Process command line for tiling.py.')
    parser.add_argument('--bed', '-b',
                        default='regions.txt', type=str, required=True, dest='regions_file',
                        help="Bed file with regions")
    parser.add_argument('--output', '-o', type=str, required=True, dest='output_file',
                        help="Output prboes file")
    parser.add_argument('--coverage', '-cov', type=str, required=True, dest='coverage_file',
                        help="Coverage file (segmented)")

    args = parser.parse_args()
    bed_file = args.regions_file
    output_file = args.output_file
    coverage_file = args.coverage_file
    merge_coverage_file(bed_file, coverage_file, output_file)


main()
