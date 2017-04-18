#!/bin/bash -e

if [ $# -lt 3 ]
then
	echo usage $0 [BAM file] [output] [whitelist]
	exit 1
fi

input=$1
output_bcf=${2}.white.bcf
output_vcf=${2}.white.vcf

whitelist=$3
#'/data/D161740/Reference/WhiteList/whitelist_20170109.txt'

ref='/data/D161740/Reference/Human/b37/human_g1k_v37.fasta'

samtools_path='/data/D161740/apps/ETC/samtools/samtools-1.3'

$samtools_path/samtools mpileup\
	-f $ref\
	-l $whitelist\
	-go $output_bcf\
	$input
	
bcftools call\
	-mO v\
	-o $output_vcf\
	$output_bcf
