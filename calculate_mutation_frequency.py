#! /usr/bin/python
import optparse,sys
import subprocess
import shutil
import glob
from utils import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input maf')
parser.add_option('-t', dest='target',help='exon target interval')
parser.add_option('-m', dest='metrics1',help='metrics file1')
parser.add_option('-M', dest='metrics2',help='metrics file2')

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
check_option(options.input, 'check options')

total_len = 0
for y in open(options.target).xreadlines():
	if y[0] != '@':
		y_1 = y.strip().split('\t')
		start = int(y_1[1])
		end = int(y_1[2])
		interval_len = end - start
		total_len += interval_len
print total_len


mutation_count = 0
for x in open(options.input).xreadlines():
	if x[0:11] != 'Hugo_Symbol':
		if x != '\n':
			mutation_count += 1

print mutation_count

frequency = (float(1000000) / float(total_len)) * float(mutation_count)

frequency_fin = '%.1f' % frequency


metrics1_bak = '%s.bak2' % options.metrics1

files1 = glob.glob(metrics1_bak)

if len(files1) != 1:
	shutil.copy2(options.metrics1, metrics1_bak)


metrics2_bak = '%s.bak2' % options.metrics2

files2 = glob.glob(metrics2_bak)

if len(files2) != 1:
	shutil.copy2(options.metrics2, metrics2_bak)

wfile1 = open(options.metrics1,'w')
a = 1
for x in open(metrics1_bak).xreadlines():
	if a == 7:
		x_1 = x.strip().split('\t')
		x_1.append('Mutation_frequency')

		wfile1.write('%s\n' % '\t'.join(x_1))
		a += 1
	elif a == 8:
		x_1 = x.strip().split('\t')
		x_1.append(frequency_fin)
		wfile1.write('%s\n' % '\t'.join(x_1))
		a += 1
	elif a < 7:
		wfile1.write(x)
		a += 1
	else:
		a += 1
		pass
wfile1.close()

wfile2 = open(options.metrics2,'w')
a = 1
for x in open(metrics2_bak).xreadlines():
	if a == 7:
		x_1 = x.strip().split('\t')
		x_1.append('Mutation_frequency')

		wfile2.write('%s\n' % '\t'.join(x_1))
		a += 1
	elif a == 8:
		x_1 = x.strip().split('\t')
		x_1.append(frequency_fin)
		wfile2.write('%s\n' % '\t'.join(x_1))
		a += 1
	elif a < 7:
		wfile2.write(x)
		a += 1
	else:
		a += 1
		pass
wfile2.close()


