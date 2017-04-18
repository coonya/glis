#!/bin/bash -e

if [ $# -lt 5 ]
then
	echo usage $0 [BAM] [target] [Output prefix] [Threads] [path of normal files]
	echo "suffix example: *T.agg.dedup.cleaned.bam --> agg.dedup.cleaned.bam"
	echo "target: SS50Mb, OPv2"
	exit 1
fi


bam=$1
target_bed=$2
output_prefix=$3
threads=$4
normal=$5

#echo Target is $target
#if [[ "$target" == 'SS50Mb']]
#then
#target_bed='/data/D161740/Reference/Target_intervals/OP_AMCv2/OP_AMCv2_baits.cnvkit.bed'

#elif [[ "$target" == 'OPv2']]
#then
#	target_bed='/data/D161740/Reference/Target_intervals/OncoPanel_v2/b37/OPv2_cnv2.bed'
#else
#	echo "Please check the target"
#	exit 1
#fi


reference='/data/D161740/Reference/Human/cnvkit/human_g1k_v37.fasta'
access_bed='/data/D161740/Reference/Human/cnvkit/access-5kb.b37.bed'
annotate_file='/data/D161740/Reference/Human/cnvkit/refFlat2.txt'
#cnvkit_PATH='/data/D161740/apps/CNV/cnvkit/cnvkit'
cnvkit_PATH='/data/D161740/apps/CNV/cnvkit/cnvkit_0.7.10'

#source /data/D161740/apps/Annotation/Oncotator/oncotator_env2/bin/activate
source /data/D161740/apps/PYTHON/venv/cnvkit/bin/activate

echo '01. running batch'
#$cnvkit_PATH/cnvkit.py batch\
cnvkit.py batch\
	$bam\
	--normal $normal\
	--targets $target_bed\
	--fasta $reference\
	--split\
	--access $access_bed\
	--output-reference reference.cnn\
	--output-dir $output_prefix\
	-p $threads $option1
#--annotate $annotate_file

for x in ${output_prefix}/*.cns
do
	sample_cns=$x
done

echo '02. call'
echo $sample_cns
echo ${sample_cns/.cns/.call.cns}

cnvkit.py call\
	$sample_cns\
	-y\
	-m clonal\
	--purity 0.5\
	-o ${sample_cns/.cns/.call.cns}

echo '04. filter out Background in cnr file'
grep -v Background ${sample_cns/.cns/.cnr} | cut -f 1,2,3,4,6 > ${sample_cns/.cns/.filtered.cnr}


echo '03. filter'
/data/D161740/src/glis/filter_cnvkit.py -i ${sample_cns/.cns/.call.cns}


deactivate
