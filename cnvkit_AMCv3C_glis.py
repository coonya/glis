#! /usr/bin/python
import optparse, os, sys, subprocess, glob
from subprocess import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='Input BAM file')
parser.add_option('-t', dest='target',help='Target bed file')
parser.add_option('-o', dest='output',help='Output dir')
parser.add_option('-T', dest='threads',help='The number of threads')
parser.add_option('-g', dest='gender',help='Gender')

(options, args) = parser.parse_args()

def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())

def check_dir(dir):
	if os.path.exists(dir):
		pass
	else:
		try:
			os.mkdir(dir)
		except:
			pass


def execute_in_virtualenv(virtualenv_name, commands):
	'''Execute Python code in a virtualenv, return its stdout and stderr.'''
	command_template = '/bin/bash -c "source {}/{  }/bin/activate && python -"'
	command = shlex.split(command_template.format(os.environ['WORKON_HOME'], virtualenv_name))
	process = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=False)
	return process.communicate(commands)

def __main__():
	check_option(options.input, "You need to use '-i' option with input vcf file\n")
	check_option(options.target, "You need to use '-t' option with input vcf file\n")
	check_option(options.output, "You need to use '-o' option with input vcf file\n")
	check_option(options.threads, "You need to use '-T' option with input vcf file\n")
	check_option(options.gender, "You need to use '-g' option with input vcf file\n")

	reference='/data/D161740/Reference/Human/cnvkit/human_g1k_v37.fasta'
	access_bed='/data/D161740/Reference/Human/cnvkit/access-5kb.b37.bed'
	annotate_file='/data/D161740/Reference/Human/cnvkit/refFlat2.txt'
	cnvkit_reference='/data/D161740/src/pipeline/ETC/amcv2/OP_AMCv2_reference.cnn'
	cnvkit_PATH='/data/D161740/apps/CNV/cnvkit/cnvkit_0.7.10'
	
	normal_bams = '/data/D161740/src/pipeline/ETC/amcv3c/*N.fastq.gz.*.bam'

	#### python virtual environment
	cmd = 'source /data/D161740/apps/PYTHON/venv/cnvkit/bin/activate'
	print 'CMD: %s' % cmd
	#subprocess.call(cmd, shell=True)

	### cnvkit batch
	cmd1 = []
	cmd1.append('/data/D161740/apps/PYTHON/venv/cnvkit/bin/python /data/D161740/apps/PYTHON/venv/cnvkit/bin/python/cnvkit.py batch')
	cmd1.append(options.input)
	cmd1.append('--normal %s' % normal_bams)
	cmd1.append('--targets %s' % options.target)
	cmd1.append('--fasta %s' % reference)
	cmd1.append('--split')
	cmd1.append('--access %s' % access_bed)
	cmd1.append('--output-reference reference.cnn')
	cmd1.append('--output-dir %s' % options.output)
	cmd1.append('-p %s' % options.threads)

	if options.gender == 'Male':
		cmd1.append('--male-reference')

	cmd = ' '.join(cmd1)
	print 'CMD: %s' % cmd
	subprocess.call(cmd, shell=True)


	cns_file = glob.glob('*.cns')

	if len(cns_file) == 1:
		if '.fastq.gz.initialAlign.merged.dedup.realign' in cns_file[0]:
			cns = cns_file[0].replace('.fastq.gz.initialAlign.merged.dedup.realign', '')
			cnr = cns_file[0].replace('.fastq.gz.initialAlign.merged.dedup.realign', '').replace('.cns', '.cnr')
			call_cns = cns.replace('.cns', '.call.cns')
			filtered_cnr = cnr.replace('.cnr', '.filtered.cnr')
	else:
		print cns_file
		sys.exit('please check the number of .cns files')
	"""	
	### CN call
	cmd2.append('cnvkit.py call')
	cmd2.append(cns)
	cmd2.append('-y')
	cmd2.append('-m clonal')
	cmd2.append('--purity 0.5')
	cmd2.append('-o %s' % call_cns)
	cmd = ' '.join(cmd2)
	print 'CMD: %s' % cmd
	subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)


	## filter out Background in .cnr file
	cmd = 'grep -v Background %s |cut -f 1,2,3,4,6 > %s' % (cnr, filtered_cnr)
	print 'CMD: %s' % cmd
	subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)


	## filter 2 copy evnets
	cmd = '/data/D161740/src/glis/filter_cnvkit.py -i %s' % call_cns
	print 'CMD: %s' % cmd
	subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)

	## deactivate python virtual environment'
	cmd = 'deactivate'
	print 'CMD: %s' % cmd
	subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)
	"""
if __name__=="__main__":__main__()


