# segmentation_before_CNV_calling
2 simple scrips which help to subdivide BED files from WES NGS results for CNV calling (works best with [ClinCNV](https://github.com/imgag/ClinCNV) )

# Motivation
Suppose we have a cohort of whole-exome sequencing results (sequenced with the same kit) and want to detect CNVs using read depth signal there. As a preparation step, we normally take a BED file with coordinates of targets and calculate coverages within the targeted regions.

## Problem 1
Targeted enrichment probes in many exome kits are of fixed length. For example, Agilent was using 120bp hybridization probes. Thus, to cover one exon of 121 base pairs length, we need at least 2 probes. Double the number of fragments - double the coverage.

## Problem 2
There are many long exons in human genome. Sometimes long variants there are not already picked up by short variant callers such as GATK (because they are too long), but also can not be detected by CNV callers since we calculate the coverage for the whole exon at once.

As a solution, we should sub-segment our BED file. However, if we simply make overlapping windows, the coverage depth signal for the overlapping windows will be heavily correlated and confuse many existing CNV callers.

# Description
This script *sub-segments targeted regions in a smart way*, imitating the probes design procedure used by Agilent, however, useful for other kits as well. We segment regions into "segment belonging to the left probe", "overlapping segment", "segment belonging to the right probe". Then we consider the coverage of the left and the right parts as "quasi-independent" from each other and divide the coverage from the overlapping segment proportionally. It solved the problems 1 and 2 with a satisfactory accuracy.

## How to run
Imagine you have a BAM file `sample.bam` and a corresponding bed file from the exome enrichment kit `kit.bed`.

One button solution (using 'one_button_run.sh' - be sure you have rights to write in this folder and run python scripts from the repository):

1) change the first line of 'one_button_run.sh' to point to list of absolute paths to BAM files you want to use (for example, as a result of 'find /path/to/your/bams/ -name "*.bam" > list_of_files.txt'). If you use CRAM files, specify reference when you run 'BedCoverage' from 'ngs-bits' (uncommend '-ref').

2) change the second line of 'one_button_run.sh' to point to your '.bed' file (3 columns - chromosome, start, end, GC or gene annotation of this file is not needed, no additional columns except 3).

3) if 'ngs-bits' is not added to 'PATH', you should change the line 'BedCoverage' to the actual '/path/to/ngs-bits/BedCoverage'. 

4) Run 'bash one_button_run.sh' and let me know if any errors appear during the run, we will try to fix them.

*As a result* you will have a folder 'prepared_coverage' which you will use later with 'ClinCNV' merging routine, and a new bed file 'segmented.name_of_your_old_bed.bed'. You should further work with this new bed file (GC-annotation, genes annotation) and do not use the original!

Manual run: at first, segment the bed file with 'probes_from_bed.py':

```
probes_from_bed.py --bed kit.bed --output supersegmented.bed --probLen 120
```

Then you calculate coverage using the `supersegmented.for_coverage.bed` file and `sample.bam`, using, for example, [ngs-bits](https://github.com/imgag/ngs-bits). After that you need to "reassemble" the coverage for `supersegmented.bed` file (use this as an input for your CNV caller instead of your original file!)

To do this you need to run `merge_segmented_coverage.py`:

```
merge_segmented_coverage.py --bed supersegmented.bed --output final_assembled_coverage.cov --coverage your_existing_coverage_file_from_the_previous_step.cov
```
