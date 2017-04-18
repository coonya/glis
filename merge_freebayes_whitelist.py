#! /usr/bin/python
import optparse, os, sys, subprocess, re, MySQLdb, glob, shutil

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input', help='freebayes output MAF')
parser.add_option('-w', dest='white', help='white list VCF')


(options, args) = parser.parse_args()
def check_dir(check_dir):
	if os.path.exists(check_dir):
		pass
	else:
		try:
			os.makedirs(check_dir)
		except:
			pass

def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())

def change_base(alt):
	base = {'A':'T', 'G':'C', 'C':'G', 'T':'A'}
	return base[alt]


def __main__():
	check_option(options.input, "You need to use '-d' option with basecall dir\n")
	check_option(options.white, "You need to use '-d' option with basecall dir\n")
	dic = {}

	out_name = options.white.replace('.vcf', '.modified.vcf')
	wfile = open(out_name, 'w')

	for x in open(options.input).xreadlines():
		if x[0:11] != 'Hugo_Symbol':
			x_1 = x.strip().split('\t')
			chr = x_1[4]
			pos = x_1[5]
			key1 = '%s_%s' % (chr, pos)
			dic[key1] = x

	for y in open(options.white).xreadlines():
		if y[0] == '#':
			wfile.write(y)
		else:
			print y
			y_1 = y.strip().split('\t')
			chr2 = y_1[0]
			pos2 = y_1[1]
			ref2 = y_1[3]
			if y_1[4] == '.':
				alt2 = change_base(ref2)
			else:
				alt2 = y_1[4]
			DP = y_1[7].split(';')[0].replace('DP=','')

			key2 = '%s_%s' % (chr2, pos2)

			if dic.has_key(key2):
				pass
			else:
				line = '%s\t%s\t.\t%s\t%s\t.\t.\t.\t%s\tGT:DP:AD\t0/0:%s:%s,0\n' % (chr2, pos2, ref2, alt2, y_1[7], DP, DP)
				wfile.write(line)

	wfile.close()





if __name__=="__main__":__main__()
