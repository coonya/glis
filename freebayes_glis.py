#! /usr/bin/python
import optparse,sys
import subprocess
from utils import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input BAM')
parser.add_option('-t', dest='target',help='target interval')
parser.add_option('-o', dest='output',help='output dir')
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


check_option(options.input, "You need to use '-i' option.")

params = read_config(options.ini)

ref = params['ref']

java_path = params['java8_path']
gatk2_path = params['gatk2_path']
vcffilter_path = params['vcffilter_path']
freebayes_path = params['freebayes_path']


output = options.input.split('/')[-1].split('.')[0]

## freebayes
cmd = '%s/freebayes -b %s -v %s/%s.freebayes.vcf -f %s -t %s -F 0.05 -C 5 -m 30 -q 20 --min-coverage 30' % (freebayes_path, options.input, options.output, output, ref, options.target)

print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout = subprocess.PIPE)

### filter variants
cmd = '%s/vcffilter -f "SAF > 0 & SAR > 0 & RPR > 1 & RPL >1" %s/%s.freebayes.vcf > %s/%s.freebayes.filtered.vcf' % (vcffilter_path, options.output, output, options.output, output)
print 'CMD: %s' % cmd
subprocess.call(cmd, shell=True, stdout = subprocess.PIPE)

#cmd = '/data/D161740/src/vcf2maf.py -i %s/%s.freebayes.filtered.vcf' % (options.output, output)
filtered_vcf = '%s/%s.freebayes.filtered.vcf' % (options.output, output)
output_maf = filtered_vcf.replace('.vcf','.maf')

line_num = 0

for k in open(filtered_vcf).xreadlines():
	if k[0] != '#':
		line_num += 1
wfile = open('%s/%s.freebayes.filtered.fixed.maf' % (options.output, output),'w')

if line_num != 0:

	cmd = '/data/D161740/src/vcf2maf.py -i %s/%s.freebayes.filtered.vcf' % (options.output, output)
	#cmd = '/data/D161740/apps/Annotation/vcf2maf/vcf2maf/vcf2maf.pl --input-vcf %s --output-maf %s --tumor-id %s' % (filtered_vcf, output_maf,
	print 'CMD: %s' % cmd
	subprocess.call(cmd, shell=True, stdout = subprocess.PIPE)


	for x in open('%s/%s.freebayes.filtered.maf' % (options.output, output)).xreadlines():
		x_1 = x.strip().split('\t')
		if x_1[0] == 'Hugo_Symbol':
			mutation_status = x_1.index('Mutation_Status')
			wfile.write(x)
		else:
			x_1[mutation_status] = 'Germline'
			x_1.append('No')
			x_1.append('No')
			wfile.write('%s\n' % '\t'.join(x_1))
else:
	wfile.write('Hugo_Symbol')

wfile.close()
