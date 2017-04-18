#! /usr/bin/python
import argparse,sys, time, os
import subprocess
from datetime import date

## Global parameters
usage = '%(prog)s [options]'
parser = argparse.ArgumentParser(prog='trimgalore.py', description='Trim adaptor sequence', usage=usage)
parser.add_argument('-f1', action='store', dest='fastq1', required=True, metavar='fastq R1', help='fastq R1 file')
parser.add_argument('-f2', action='store', dest='fastq2', required=True, metavar='fastq R2', help='fastq R2 file')
parser.add_argument('-o', action='store', dest='output', required=True, metavar='output', help='output dir')

#parser.add_argument('-a', action='store', dest='barcode', required=True, metavar='1', help='barcode number')
args = parser.parse_args()


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
	

start_time = time.time()

output = args.output

if not os.path.exists(output):
	os.makedirs(output)

trimgalore_path = '/data/D161740/apps/ETC/trimgalore/trim_galore_zip'
cutadapt_path = '/data/D161740/apps/ETC/cutadapt/cutadapt-1.9/bin/cutadapt'

#adaptor1 = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC%sATCTCGTATGCCGTCTTCTGCTTG' % barcode(args.barcode)
adaptor1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
adaptor2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'


#cmd = '%s/trim_galore --path_to_cutadapt %s --paired --gzip -q 1 --suppress_warn --stringency 3 -length 25 -o %s -a %s -a2 %s %s %s' % (trimgalore_path, cutadapt_path, output, barcode(args.barcode), adaptor2, args.fastq1, args.fastq2)

### command line

cmd = []
cmd.append('%s/trim_galore' % trimgalore_path)
cmd.append('--path_to_cutadapt %s' % cutadapt_path)
cmd.append('--paired')
cmd.append('--gzip')
cmd.append('-q 1')
cmd.append('--suppress_warn')
cmd.append('--stringency 3')
cmd.append('--length 25')
cmd.append('-o %s' % output)
cmd.append('-a %s' % adaptor1)
cmd.append('-a2 %s' % adaptor2)
cmd.append('-q 5')
cmd.append('%s' % args.fastq1)
cmd.append('%s' % args.fastq2)


final_cmd = '\t'.join(cmd)
print 'CMD: %s\n' % final_cmd


#subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)
subprocess.call(final_cmd, shell=True)

end_time = time.time()


#print 'elepsed time was %g.' % (
