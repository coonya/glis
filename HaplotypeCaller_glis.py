#! /usr/bin/python
import optparse,sys
import subprocess
from utils import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input BAM')
parser.add_option('-t', dest='target',help='target interval')
parser.add_option('-o', dest='output',help='output')
parser.add_option('-c', dest='ini',help='ini file')

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

params = read_config(options.ini)

ref = params['ref']

java_path = params['java8_path']
gatk2_path = params['gatk2_path']


## running Haplotype Caller
cmd = '%s/java -jar %s/GenomeAnalysisTK.jar -T HaplotypeCaller -R %s -I %s -L %s -o %s' % (java_path, gatk2_path, ref, options.input, options.target, options.output)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout = subprocess.PIPE)


output_maf = options.output.replace('.vcf', '.maf')
cmd = 'vcf2maf.pl --input-vcf %s --output-maf %s --tumor-id %s' % (wname, output_maf, sample)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)



wname = output_maf.replace('.maf', '.fixed.maf')
wfile = open(wname, 'w')

for x in open(output_maf).xreadlines():
	x_1 = x.replace('\n','').split('\t')
	if x_1[0] == 'Hugo_Symbol':
		mutation_status = x_1.index('Mutation_Status')
		wfile.write(x)
	
	else:
		x_1[mutation_status] = 'Germline'
		wfile.write('%s\n' % '\t'.join(x_1))
wfile.close()
