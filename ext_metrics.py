#! /usr/bin/python
import optparse,sys
import subprocess
from utils import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input')
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

wfile = open(options.output, 'w')
a = 1
for x in open(options.input).xreadlines():
	if a == 7:
		x_1 = x.strip().split('\t')
		bait_set = x_1.index('BAIT_SET')
		bait_territory = x_1.index('BAIT_TERRITORY')
		target_territory = x_1.index('TARGET_TERRITORY')
		total_reads = x_1.index('TOTAL_READS')
		pf_reads = x_1.index('PF_READS')
		pf_unique_reads = x_1.index('PF_UNIQUE_READS')
		pct_pf_reads = x_1.index('PCT_PF_READS')
		pct_pf_uq_reads = x_1.index('PCT_PF_UQ_READS')
		pct_selected_bases = x_1.index('PCT_SELECTED_BASES')
		pct_off_bait = x_1.index('PCT_OFF_BAIT')
		mean_target_coverage = x_1.index('MEAN_TARGET_COVERAGE')
		median_target_coverage = x_1.index('MEDIAN_TARGET_COVERAGE')
		zero_cvg_targets_pct = x_1.index('ZERO_CVG_TARGETS_PCT')
		pct_target_bases_30x = x_1.index('PCT_TARGET_BASES_30X')
		pct_target_bases_100x = x_1.index('PCT_TARGET_BASES_100X')

		line = []

		line.append(x_1[bait_set])
		line.append(x_1[bait_territory])
		line.append(x_1[target_territory])
		line.append(x_1[total_reads])
		line.append(x_1[pf_reads])
		line.append(x_1[pf_unique_reads])
		line.append(x_1[pct_pf_reads])
		line.append(x_1[pct_pf_uq_reads])
		line.append(x_1[pct_selected_bases])
		line.append(x_1[pct_off_bait])
		line.append(x_1[mean_target_coverage])
		line.append(x_1[median_target_coverage])
		line.append(x_1[zero_cvg_targets_pct])
		line.append(x_1[pct_target_bases_30x])
		line.append(x_1[pct_target_bases_100x])

		wfile.write('%s\n' % '\t'.join(line))
		a += 1
	elif a == 8:
		x_1 = x.strip().split('\t')
		line = []

		line.append(x_1[bait_set])
		line.append(x_1[bait_territory])
		line.append(x_1[target_territory])
		line.append(x_1[total_reads])
		line.append(x_1[pf_reads])
		line.append(x_1[pf_unique_reads])
		line.append(x_1[pct_pf_reads])
		line.append(x_1[pct_pf_uq_reads])
		line.append(x_1[pct_selected_bases])
		line.append(x_1[pct_off_bait])
		line.append(x_1[mean_target_coverage])
		line.append(x_1[median_target_coverage])
		line.append(x_1[zero_cvg_targets_pct])
		line.append(x_1[pct_target_bases_30x])
		line.append(x_1[pct_target_bases_100x])

		wfile.write('%s\n' % '\t'.join(line))
		a += 1
	elif a < 7:
		wfile.write(x)
		a += 1
	else:
		a += 1
		pass
wfile.close()



