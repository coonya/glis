#!/bin/bash -e

if [ $# -lt 4 ]
then
	echo usage $0 [BAM] [REF] [target_interval] [output_dir]
	exit 1
fi

vcffilter_path=/data/D161740/apps/ETC/vcflib/vcflib/bin
freebayes_path=/data/D161740/apps/Cancer/freebayes/bin
ref=$2
target_interval=$3
output_dir=$4

tumor_bam=$1

output=$(echo $tumor_bam|rev|cut -d '/' -f1|rev|cut -d '.' -f1)
echo $output

## freebayes
$freebayes_path/freebayes\
	-b $tumor_bam\
	-v $output_dir/${output}.freebayes.vcf\
	-f $ref\
	-t $target_interval\
	-F 0.05\
	-C 5\
	-m 30\
	-q 20\
	--min-coverage 30

### filter variants
$vcffilter_path/vcffilter\
	-f "SAF > 0 & SAR > 0 & RPR > 1 & RPL >1"\
	$output_dir/${output}.freebayes.vcf\
	> $output_dir/${output}.freebayes.filtered.vcf

/data/D161740/src/vcf2maf.py -i $output_dir/${output}.freebayes.filtered.vcf


