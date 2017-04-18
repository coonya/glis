#! /usr/bin/python
import optparse,sys
import subprocess
import shutil
import glob
from utils import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input')
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
input_bak = '%s.bak' % options.input
files = glob.glob(input_bak)

if len(files) != 1:
	shutil.copy2(options.input, input_bak)


wfile = open(options.input, 'w')
a = 1
for x in open(input_bak).xreadlines():
	print a, x
	if a == 7:
		x_1 = x.strip().split('\t')
		category = x_1.index('CATEGORY')
		total_reads = x_1.index('TOTAL_READS')
		pf_reads = x_1.index('PF_READS')
		pct_pf_reads = x_1.index('PCT_PF_READS')
		pf_reads_aligned = x_1.index('PF_READS_ALIGNED')
		pct_pf_reads_aligned = x_1.index('PCT_PF_READS_ALIGNED')
		pf_aligned_bases = x_1.index('PF_ALIGNED_BASES')
		pf_mismatch_rate = x_1.index('PF_MISMATCH_RATE')
		pf_indel_rate = x_1.index('PF_INDEL_RATE')
		mean_read_length = x_1.index('MEAN_READ_LENGTH')
		reads_algined_in_pairs = x_1.index('READS_ALIGNED_IN_PAIRS')
		pct_reads_aligned_in_pairs = x_1.index('PCT_READS_ALIGNED_IN_PAIRS')
		pct_chimeras = x_1.index('PCT_CHIMERAS')
		pct_adapter = x_1.index('PCT_ADAPTER')

		line = []

		line.append(x_1[category])
		line.append(x_1[total_reads])
		line.append(x_1[pf_reads])
		line.append(x_1[pct_pf_reads])
		line.append(x_1[pf_reads_aligned])
		line.append(x_1[pct_pf_reads_aligned])
		line.append(x_1[pf_aligned_bases])
		line.append(x_1[pf_mismatch_rate])
		line.append(x_1[pf_indel_rate])
		line.append(x_1[mean_read_length])
		line.append(x_1[reads_algined_in_pairs])
		line.append(x_1[pct_reads_aligned_in_pairs])
		line.append(x_1[pct_chimeras])
		line.append(x_1[pct_adapter])

		wfile.write('%s\n' % '\t'.join(line))
		a += 1
	elif a >=8 and a <= 10:
		x_1 = x.strip().split('\t')
		line = []

		line.append(x_1[category])
		line.append(x_1[total_reads])
		line.append(x_1[pf_reads])
		line.append(x_1[pct_pf_reads])
		line.append(x_1[pf_reads_aligned])
		line.append(x_1[pct_pf_reads_aligned])
		line.append(x_1[pf_aligned_bases])
		line.append(x_1[pf_mismatch_rate])
		line.append(x_1[pf_indel_rate])
		line.append(x_1[mean_read_length])
		line.append(x_1[reads_algined_in_pairs])
		line.append(x_1[pct_reads_aligned_in_pairs])
		line.append(x_1[pct_chimeras])
		line.append(x_1[pct_adapter])

		wfile.write('%s\n' % '\t'.join(line))
		a += 1
	elif a < 7:
		wfile.write(x)
		a += 1
	else:
		a += 1
		pass
wfile.close()



