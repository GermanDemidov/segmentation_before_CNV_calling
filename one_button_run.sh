list_of_files="/path/to/list_of_files.txt"
bed_file="/path/to/bed/file.bed"

final_output_dir="./prepared_coverage/"
tmp_output_dir="./tmp/"


echo "Creating folders. If the scripts fails here, you don't have permissions to create these folders here. Do chmod on this folder."
mkdir $tmp_output_dir
mkdir $final_output_dir

./probes_from_bed.py --bed $bed_file --output ./tmp/supersegmented.bed --probLen 120
sort -k1,1V -k2,2n -k3,3n  ./tmp/supersegmented.for_coverage.bed > ./tmp/sorted.supersegmented.for_coverage.bed

s=${bed_file##*/}
bed_file_supersegmented="segmented."$s

echo "Your future bed file for CNV calling is "$bed_file_supersegmented

while IFS= read -r line
do
    echo "Calculating coverage for file: $line"
	s=${line##*/}
	name="${s%.*}"
	BedCoverage -bam $line -in ./tmp/sorted.supersegmented.for_coverage.bed -out $tmp_output_dir"/"$name".cov" -decimals 4 -min_mapq 2 -threads 4 # -ref PATH_TO_REFERENCE_IF_YOU_USE_CRAM!
    ./merge_segmented_coverage.py --bed ./tmp/supersegmented.bed --output $final_output_dir"/"$name".cov" --coverage $tmp_output_dir"/"$name".cov"
done < $list_of_files

cp ./tmp/supersegmented.bed ./$bed_file_supersegmented
rm -r ./tmp/

