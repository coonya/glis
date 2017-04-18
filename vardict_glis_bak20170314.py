#! /usr/bin/python
import optparse,sys
import subprocess

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input BAM')
parser.add_option('-t', dest='threads',help='The number of threads')
parser.add_option('-L', dest='target',help='target bed')
parser.add_option('-o', dest='output',help='output')
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

check_option(options.input, "You need to use '-i' option.")



ref="/data/D161740/Reference/Human/b37/human_g1k_v37.fasta"
tumor = options.input
tumor_name = tumor.split('.')[0]
target = options.target
threads = options.threads

vardict_path="/data/D161740/apps/Cancer/VarDict/VarDict-1.4.10/bin"
vardict_path2="/data/D161740/apps/Cancer/VarDict/VarDict"


AF_THR="0.05"

vardict_out = options.output

#vardict_out2 = options.output.replace('.vcf', '.tmp.vcf')
## running Vardict
#cmd = '%s/VarDict -C -G %s -t -f %s -N %s -b %s -c 1 -S 2 -E 3 -g 5 -th %s %s | %s/teststrandbias.R | %s/var2vcf_valid.pl -N %s -f %s -Q 20 -d 30 -v 3 -c 10 > %s' % (vardict_path, ref, AF_THR, tumor_name, tumor, threads, target, vardict_path2, vardict_path2, tumor_name, AF_THR, vardict_out)
#cmd = '%s/VarDict -r 3 -C -G %s -t -f %s -N %s -b %s -c 1 -S 2 -E 3 -g 5 -th %s %s | %s/teststrandbias.R | %s/var2vcf_valid.pl -N %s -f %s -Q 20 -d 30 -v 3 -c 10 > %s' % (vardict_path, ref, AF_THR, tumor_name, tumor, threads, target, vardict_path2, vardict_path2, tumor_name, AF_THR, vardict_out)
cmd = '%s/VarDict -C -G %s -t -f %s -N %s -b %s -c 1 -S 2 -E 3 -g 5 -th %s %s | %s/teststrandbias.R | %s/var2vcf_valid.pl -N %s -f %s -Q 20 -d 30 -v 3 -c 10 > %s' % (vardict_path, ref, AF_THR, tumor_name, tumor, threads, target, vardict_path2, vardict_path2, tumor_name, AF_THR, vardict_out)
#cmd = '%s/VarDict -G %s -f %s -N %s -b %s -c 1 -S 2 -E 3 -g 5 -th %s %s | %s/teststrandbias.R | %s/var2vcf_valid.pl -N %s -f %s > %s' % (vardict_path, ref, AF_THR, tumor_name, tumor, threads, target, vardict_path2, vardict_path2, tumor_name, AF_THR, vardict_out)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout = subprocess.PIPE)

#cmd = 'vcf-sort -c %s > %s' % (vardict_out2, vardict_out)
#print 'CMD: %s' % cmd
#subprocess.call(cmd, shell=True, stdout = subprocess.PIPE)


## filtering
filtered_vcf = vardict_out.replace('.vcf', '.filtered.vcf')
wfile = open(filtered_vcf, 'w')

for x in open(vardict_out).xreadlines():
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

cmd = '/data/D161740/src/filter_SNV.py -i %s' % anno_maf
#cmd = '/data/D161740/src/glis/filter_variants.py -i %s -d %s -D %s -p %s -P %s -k %s -K %s -r %s -b %s' % (anno_maf, options.dbsnp_snv, options.dbsnp_indel, options.pon_snv, options.pon_indel, options.KRG_common, options.KRG_indel, options.KRG_rare, options.black)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)
