#! /usr/bin/python
import optparse, os, sys, subprocess, re, MySQLdb, glob, shutil

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-s', dest='input', help='vardict output MAF')
parser.add_option('-w', dest='white', help='white list maf')
parser.add_option('-o', dest='output', help='output')


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

def __main__():
	check_option(options.input, "You need to use '-s' option with basecall dir\n")
	check_option(options.white, "You need to use '-w' option with basecall dir\n")
	dic = {}

	wfile = open(options.output, 'w')

	for x in open(options.input).xreadlines():
		x_1 = x.strip().split('\t')
		key1 = '%s_%s' % (x_1[4], x_1[5])
		dic[key1] = x
		wfile.write(x)


	for y in open(options.white).xreadlines():
		y_1 = y.strip().split('\t')
		key2 = '%s_%s' % (y_1[4], y_1[5])

		if len(dic) == 1:
			wfile.write(dic['%s_%s' % (Chromosome, Start_Position)])
		else:
			if dic.has_key(key2):
				pass
			else:
				wfile.write(y)

	wfile.close()





if __name__=="__main__":__main__()
