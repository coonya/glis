#!/bin/bash -e
if [ $# -lt 3 ]
then
	echo usage: $0 [BAM] [target] [output]
	exit 1
fi

### variables

ref='/data/D161740/Reference/Human/b37/human_g1k_v37.fasta'
gatk_path='/data/D161740/apps/ETC/GATK/GenomeAnalysis-3.5.0'

bam=$1
target=$2
output=$3



/usr/bin/java\
	-jar $gatk_path/GenomeAnalysisTK.jar\
	-T HaplotypeCaller\
	-R $ref\
	-I $bam\
	-L $target\
	-o $output
