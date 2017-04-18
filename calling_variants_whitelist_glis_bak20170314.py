#! /usr/bin/python
import optparse,sys
import subprocess

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input BAM')
parser.add_option('-w', dest='white',help='white list')
parser.add_option('-o', dest='output',help='output')

(options, args) = parser.parse_args()


##### definitions
def logo():

	logo = """
###############################################################################
### @@@@@@@ ###################################################################
## @@@    @@ ###     #####     #####     ##############     ###################
## @@ ###    ## @@@@@ ### @@@@@ ### @@@@@ ## @@ # @@ # @@@@@@ #################
## @@ ######## @@ # @@ # @@ # @@ # @@   @@ ## @@ @@ # @@    @@ ################
## @@ ###    # @@ # @@ # @@ # @@ # @@ # @@ ### @@@ # @@ #### @@ ###############
## @@@    @@ # @@ # @@ # @@ # @@ # @@ # @@ ### @@ ### @@    @@@@ ##############
### @@@@@@@ ### @@@@@ ### @@@@@ ## @@ # @@ ## @@ ##### @@@@@ # @@ #############
###############################################################################
######################################### OncoPanelV2 Analysis Pipeline v2.0 ##
###############################################################################
"""
	print logo
	


def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())

def change(base):
	dic = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	return dic[base]

def change_line(line):
	x_1 = line.strip().split('\t')
	ref = x_1[3]
	x_1[4] = change(ref)

	x_2 = x_1[7].split(';')
	DP = x_2[0].replace('DP=','')

	x_1[8] = 'GT:DP:AD'
	x_1[9] = '0/1:%s:%s,0' % (DP, DP)

	fixed_line = '\t'.join(x_1)

	return fixed_line



check_option(options.input, "You need to use '-i' option.")



ref = "/data/D161740/Reference/Human/b37/human_g1k_v37.fasta"

samtools_path = '/data/D161740/apps/ETC/samtools/samtools-1.3'


output = options.output.replace('.vcf', '.bcf')

white_name = options.white.replace('.txt', '.input.txt')
white_input = open(white_name, 'w')
white_dic = {}

for l in open(options.white).xreadlines():
	l_1 = l.strip().split('\t')
	white_input.write('%s\t%s\n' % (l_1[0], l_1[1]))
	white_key = '%s_%s_%s' % (l_1[0], l_1[1], l_1[2])
	white_dic[white_key] = '1'

white_input.close()


## running Vardict
cmd = '%s/samtools mpileup -f %s -l %s -go %s %s' % (samtools_path, ref, white_name, output, options.input)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout = subprocess.PIPE)


## filtering
### annotation with VEP
output_vcf = options.output

cmd = 'bcftools call -mO v -o %s %s' % (output_vcf, output)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)

wname = output_vcf.replace('.vcf','.modified.vcf')
wfile = open(wname, 'w')

for x in open(output_vcf).xreadlines():
	if x[0] == '#':
		wfile.write(x)
		if x[0:6] == '#CHROM':
			sample = x.strip().split('\t')[9]
	else:
		x_1 = x.strip().split('\t')
		vcf_key = '%s_%s_%s' % (x_1[0], x_1[1], x_1[3])
		if white_dic.has_key(vcf_key):
			wfile.write('%s\n' % change_line(x))

wfile.close()

output_maf = wname.replace('.vcf', '.maf')
cmd = 'vcf2maf.pl --input-vcf %s --output-maf %s --tumor-id %s' % (wname, output_maf, sample)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)

wname2 = output_maf.replace('.maf', '.fixed.maf')
wfile2 = open(wname2, 'w')

for y in open(output_maf).xreadlines():
	y_1 = y.replace('\n','').split('\t')
	if y[0:11] == 'Hugo_Symbol':
		ref_pos = y_1.index('Reference_Allele')
		alt_pos = y_1.index('Tumor_Seq_Allele2')
		HGVSc_pos = y_1.index('HGVSc')
		HGVSp_pos = y_1.index('HGVSp')
		HGVSp_short_pos = y_1.index('HGVSp_Short')
		all_effect_pos = y_1.index('all_effects')
		allele_pos = y_1.index('Allele')
		consequence_pos = y_1.index('Consequence')
		aminoacid_pos = y_1.index('Amino_acids')
		codon_pos = y_1.index('Codons')
		id_pos = y_1.index('Existing_variation')
		sift_pos = y_1.index('SIFT')
		polyphen_pos = y_1.index('PolyPhen')

	else:
		y_1[alt_pos] = y_1[ref_pos]
		y_1[HGVSc_pos] = '%s%s' % (y_1[HGVSc_pos][:-1], y_1[ref_pos])
		y_1[HGVSp_pos] = '%s%s' % (y_1[HGVSp_pos][:-3], y_1[HGVSp_pos][2:5])
		y_1[HGVSp_short_pos] = '%s%s' % (y_1[HGVSp_short_pos][:-1], y_1[HGVSp_pos][2])
		y_1[all_effect_pos] = '-'
		y_1[consequence_pos] = '-'
		y_1[allele_pos] = y_1[ref_pos]
		y_1[aminoacid_pos] = y_1[aminoacid_pos].replace(y_1[aminoacid_pos][-1], y_1[aminoacid_pos][0])
		y_1[codon_pos] = y_1[codon_pos].replace(y_1[codon_pos][4:], y_1[codon_pos][0:3])
		y_1[id_pos] = '-'
		y_1[sift_pos] = '-'
		y_1[polyphen_pos] = '-'
		y_1.append('No')
		y_1.append('No')
		y_1.append('No')
		y_1.append('No')
		y_1.append('No')
		y_1.append('No')
		y_1.append('Yes')

		line = '\t'.join(y_1)

		wfile2.write('%s\n' % line)
wfile2.close()
