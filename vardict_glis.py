#! /usr/bin/python
import optparse,sys
import subprocess
from utils import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input BAM')
parser.add_option('-t', dest='threads',help='The number of threads')
parser.add_option('-L', dest='target',help='target bed')
parser.add_option('-o', dest='output',help='output')
parser.add_option('-c', dest='ini',help='ini file')




"""
parser.add_option('-d', dest='dbsnp_snv',help='dbsnp of SNV')
parser.add_option('-D', dest='dbsnp_indel',help='dbsnp of INDEL')
parser.add_option('-p', dest='pon_snv',help='PON of SNV')
parser.add_option('-P', dest='pon_indel',help='PON of INDEL')
parser.add_option('-k', dest='KRG_common',help='KRG common')
parser.add_option('-K', dest='KRG_indel',help='KRG indel')
parser.add_option('-r', dest='KRG_rare',help='KRG rare')
parser.add_option('-b', dest='black',help='black list')
"""



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

def exonic(type):
	exon_dic = {'Missense_Mutation':'', 'Nonsense_Mutation':'', 'Splice_Site':'', 'In_Frame_Del':'','Frame_Shift_Del':'','Frame_Shift_Ins':'', 'Nonstop_Mutation':'', 'In_Frame_Ins':''}
	if exon_dic.has_key(type):
		val = 'Yes'
	else:
		val = 'No'
	return val

check_option(options.input, "You need to use '-i' option.")
check_option(options.ini, "You need to use '-c' option.")

params = read_config(options.ini)


ref = params['ref']
tumor = options.input
tumor_name = tumor.split('.')[0]
target = options.target
threads = options.threads

vardict_path=params['vardict_path']
vardict_path2=params['vardict_path2']

script_bin = '/data/D161740/src/glis'

AF_THR=params['af_thr']

vardict_out = options.output

vardict_out2 = options.output.replace('.vcf', '.tmp.vcf')
## running Vardict
#cmd = '%s/VarDict -C -G %s -t -f %s -N %s -b %s -c 1 -S 2 -E 3 -g 5 -th %s %s | %s/teststrandbias.R | %s/var2vcf_valid.pl -N %s -f %s -Q 20 -d 30 -v 3 -c 10 > %s' % (vardict_path, ref, AF_THR, tumor_name, tumor, threads, target, vardict_path2, vardict_path2, tumor_name, AF_THR, vardict_out)
#cmd = '%s/VarDict -r 3 -C -G %s -t -f %s -N %s -b %s -c 1 -S 2 -E 3 -g 5 -th %s %s | %s/teststrandbias.R | %s/var2vcf_valid.pl -N %s -f %s -Q 20 -d 30 -v 3 -c 10 > %s' % (vardict_path, ref, AF_THR, tumor_name, tumor, threads, target, vardict_path2, vardict_path2, tumor_name, AF_THR, vardict_out)
#cmd = '%s/VarDict -C -G %s -t -f %s -N %s -b %s -c 1 -S 2 -E 3 -g 5 -th %s %s | %s/teststrandbias.R | %s/var2vcf_valid.pl -N %s -f %s -Q 20 -d 30 -v 3 -c 10 > %s' % (vardict_path, ref, AF_THR, tumor_name, tumor, threads, target, vardict_path2, vardict_path2, tumor_name, AF_THR, vardict_out2)
cmd = '%s/VarDict -C -G %s -t -N %s -b %s -c 1 -S 2 -E 3 -g 5 -th %s %s | %s/teststrandbias.R | %s/var2vcf_valid.pl -N %s -f %s -Q 20 -d 30 -v 3 -c 10 > %s' % (vardict_path, ref, tumor_name, tumor, threads, target, vardict_path2, vardict_path2, tumor_name, AF_THR, vardict_out2)
#cmd = '%s/VarDict -G %s -f %s -N %s -b %s -c 1 -S 2 -E 3 -g 5 -th %s %s | %s/teststrandbias.R | %s/var2vcf_valid.pl -N %s -f %s > %s' % (vardict_path, ref, AF_THR, tumor_name, tumor, threads, target, vardict_path2, vardict_path2, tumor_name, AF_THR, vardict_out)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout = subprocess.PIPE)

cmd = 'cat %s | vcf-sort -c > %s' % (vardict_out2, vardict_out)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout = subprocess.PIPE)


## filtering

cmd = '%s/vardict_filter.py -i %s -c %s' % (script_bin, vardict_out, options.ini)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout = subprocess.PIPE)

anno_vcf = vardict_out.replace('.vcf', '.anno.vcf')
filtered_vcf = vardict_out.replace('.vcf', '.anno.filtered.vcf')
wfile = open(filtered_vcf, 'w')

for x in open(anno_vcf).xreadlines():
	if x[0] == '#':
		wfile.write(x)
	else:
		x_1 = x.strip().split('\t')
		if x_1[6] == 'PASS':
			wfile.write(x)
wfile.close()

### annotation with VEP
anno_maf = filtered_vcf.replace('.vcf', '.maf')

cmd = 'vcf2maf.pl --input-vcf %s --output-maf %s --tumor-id %s' % (filtered_vcf, anno_maf, tumor_name)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)


#cmd = '/data/D161740/src/glis/add_annotation_glis.py -i %s' % anno_maf
#print 'CMD: %s' % cmd
#subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)


add_maf1 = anno_maf.replace('.maf', '.addfilter.maf')
add_maf2 = anno_maf.replace('.maf', '.addfilter.exon.maf')
wfile1 = open(add_maf1, 'w')
wfile2 = open(add_maf2, 'w')

for y in open(anno_maf).xreadlines():
	y_1 = y.replace('\n','').split('\t')
	if y_1[0] == 'Hugo_Symbol':
		tumor_id = y_1.index('Tumor_Sample_Barcode')
		AF_EAS = y_1.index('ExAC_AF_EAS')
		var_type = y_1.index('Variant_Classification')
		mutation_status= y_1.index('Mutation_Status')
		y_1.append('ExAC_filter')
		y_1.append('WhiteList')
		wfile1.write('%s\n' % '\t'.join(y_1))
		wfile2.write('%s\n' % '\t'.join(y_1))

	else:
		y_1[tumor_id] = y_1[tumor_id].split('/')[-1]
		y_1[mutation_status] = 'Somatic'
		if exonic(y_1[var_type]) == 'Yes':
			if y_1[AF_EAS] != '':
				if float(y_1[AF_EAS]) >= 0.01:
					y_1.append('Yes')
					y_1.append('No')
				else:
					y_1.append('No')
					y_1.append('No')
					wfile2.write('%s\n' % '\t'.join(y_1))

				wfile1.write('%s\n' % '\t'.join(y_1))
			else:
				y_1.append('No')
				y_1.append('No')
				wfile1.write('%s\n' % '\t'.join(y_1))
				wfile2.write('%s\n' % '\t'.join(y_1))

wfile1.close()
wfile2.close()
