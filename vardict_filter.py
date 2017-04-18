#! /usr/bin/python
import optparse, os, sys, subprocess
from subprocess import *
from utils import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input vcf file for filtering')
parser.add_option('-c', dest='ini',help='ini file')

#parser.add_option('-d', dest='dbsnp',help='dbsnp.vcf.gz')
#parser.add_option('-p', dest='pon',help='PON.vcf.gz')
#parser.add_option('-g', dest='gnomad',help='gnomad.vcf.gz')
#parser.add_option('-k', dest='krg',help='KRG.vxf.gz')

#parser.add_option('-c', dest='caller',help='somatic mutation caller e.g.) mutect, somaticindelocator')

(options, args) = parser.parse_args()

def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())

def db_com(query, dic):
	if dic.has_key(query):
		val = 'Yes'
	else:
		val = 'No'
	return val

def exac(query):
	if query == '':
		val = 'No'
	elif float(query) >= 0.01:
		val = 'Yes'
	elif float(query) < 0.01:
		val = 'No'
	else:
		print query
	return val


def read_db(db, chr, pos, ref, alt):
	cmd = 'tabix %s %s:%s-%s' % (db, chr, pos, pos)

	pipe = Popen(cmd, shell=True, stdout=PIPE)

	con = pipe.stdout.read()
	
	if con == '':
		val = 'No'
	else:
		con1 = con.split('\n')
		if len(con1) == 2:
			con2 = con.split('\t')
		
			if ',' in con2[4]:
				for y in con2[4].split(','):
	
					if con2[3] == ref and con2[4] == alt:
						val = 'Yes'
					else:
						val = 'No'
			else:
				if con2[3] == ref and con2[4] == alt:
					val = 'Yes'
				else:
					val = 'No'
		else:
			dic = {}

			for z in con1:
				con2 = con.split('\t')
				key = '%s_%s' % (chr, pos)
		
				if ',' in con2[4]:
					for y in con2[4].split(','):
	
						if con2[3] == ref and con2[4] == alt:
							dic[key] = z
				else:
					if con2[3] == ref and con2[4] == alt:
						dic[key] = z
			if dic.has_key(key):
				val = 'Yes'
			else:
				val = 'No'
	
	return val

def __main__():
	check_option(options.input, "You need to use '-i' option with input vcf file\n")
	check_option(options.ini, "You need to use '-c' option with ini file\n")

	input = options.input

	params = read_config(options.ini)

	dbsnp = params['dbsnp_common']
	pon = params['pon']
	krg = params['krg']
	gnomad = params['gnomad']
	black_list = params['black_list']



	wname = input.replace('.vcf', '.anno.vcf') 
	wfile = open(wname,'w') 
	chr_dic = {}

	for x in open(input).xreadlines(): 
		if x[0] == '#':
			wfile.write(x)
		else:
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[3]
			var = x_1[4]
			anno = x_1[7].split(';')
			type = anno[1].replace('TYPE=','')


			if chr_dic.has_key(chr):
				pass
			else:
				chr_dic[chr] = ''
				print 'Parsing chromosome%s .....' % chr

			#####
			dbsnp_com = read_db(dbsnp, chr, pos, ref, var)
			anno.append('DBSNP_COMMON=%s' % dbsnp_com)
			pon_com = read_db(pon, chr, pos, ref, var)
			anno.append('PON=%s' % pon_com)
			krg_com = read_db(krg, chr, pos, ref, var)
			anno.append('KRG=%s' % krg_com)
			gnomad_com = read_db(gnomad, chr, pos, ref, var)
			anno.append('gnomAD=%s' % gnomad_com)
			black_list_com = read_db(black_list, chr, pos, ref, var)
			anno.append('BLACKLIST=%s' % black_list_com)
			info = ';'.join(anno)
			x_1[7] = info
			if dbsnp_com == 'Yes' or pon_com == 'Yes' or krg_com == 'Yes' or gnomad_com == 'Yes' or black_list_com == 'Yes':
				x_1[6] = 'REJECT'
			else:
				if x_1[6] != 'PASS':
					x_1[6] = 'REJECT'
				else:
					x_1[6] = 'PASS'
			wfile.write('%s\n' % '\t'.join(x_1))

				
	wfile.close()
	print 'Completed!!!'

if __name__=="__main__":__main__()
