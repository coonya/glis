#! /usr/bin/python
import optparse, os, sys, subprocess
from subprocess import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input VCF file for filtering')
#parser.add_option('-c', dest='caller',help='somatic mutation caller e.g.) mutect, somaticindelocator')

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

def checkorder(line, mode1):
	x_1 = line.strip().split('\t')
	if mode1 == 'paired':
		if x_1[9][-1] == 'N' and x_1[10][-1] =='T':
			normal_col = 9
			tumor_col = 10
			mode = 'paired'
		elif x_1[9][-1] == 'T' and x_1[10][-1] =='N':
			normal_col = 10
			tumor_col = 9
			mode = 'paired'
		elif x_1[9].upper() == 'NONE':
			normal_col = 'NA'
			tumor_col = 10
			mode = 'unpaired'
		elif x_1[10].upper() == 'NONE':
			normal_col = 'NA'
			tumor_col = 9
			mode = 'unpaired'
	elif mode1 == 'unpaired':
		tumor_col = 9
		normal_col = 'NA'
		mode = 'unpaired'

	return normal_col, tumor_col, mode
	

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

def __main__():
	check_option(options.input, "You need to use '-i' option with input vcf file\n")
	input_maf = options.input

	#dbsnp_file = '/data/D161740/Reference/dbsnp/dbsnp_141_common_b37.snv.vcf'
	dbsnp_file = '/data/D161740/Reference/dbsnp/v149/dbsnp_149_common_b37.snv.vcf'
	PON_file = '/data/D161740/Reference/PON/mutect_PON_20151116_1000samples_over10.splited.vcf'
	KRG_common_file = '/data/D161740/Reference/KRG/KRG1100/KRG1100_common.txt'
	KRG_rare_file = '/data/D161740/Reference/KRG/KRG1100/KRG1100_rare.txt'

	dbsnp_dic ={}
	PON_dic = {}
	KRG_common_dic = {}
	KRG_rare_dic = {}

	for x in open(dbsnp_file):
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[3]
			alt = x_1[4]
			key = '%s_%s_%s_%s' % (chr, pos, ref, alt)
			dbsnp_dic[key] = ''
	print 'Complete reading dbsnp file.'
	
	for x in open(PON_file):
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[3]
			alt = x_1[4]
			key = '%s_%s_%s_%s' % (chr, pos, ref, alt)
			PON_dic[key] = ''

	print 'Complete reading PON file.'
	for x in open(KRG_common_file):
		if x[0:3] != 'chr':
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[4]
			alt = x_1[5]
			key = '%s_%s_%s_%s' % (chr, pos, ref, alt)
			KRG_common_dic[key] = ''

	print 'Complete reading KRG common file.'
	for x in open(KRG_rare_file):
		if x[0:3] != 'chr':
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[4]
			alt = x_1[5]
			key = '%s_%s_%s_%s' % (chr, pos, ref, alt)
			KRG_rare_dic[key] = ''
	
	print 'Complete reading KRG rare file.'
	wname = input_maf.replace('.maf', '.filter.maf') 
	wfile = open(wname,'w') 

	wname2 = input_maf.replace('.maf', '.filter.exon.maf')
	wfile2 = open(wname2,'w')

	exon_dic = {'Missense_Mutation':'', 'Nonsense_Mutation':'', 'Splice_Site':'', 'In_Frame_Del':'','Frame_Shift_Del':'','Frame_Shift_Ins':'', 'Nonstop_Mutation':'', 'In_Frame_Ins':''}
	
	for x in open(input_maf).xreadlines(): 
		if x[0:11] == 'Hugo_Symbol':
			x_1 = x.replace('\n','').split('\t')
			ExAC_AF_EAS_pos = x_1.index('ExAC_AF_EAS')

			x_1.append('Common_dbsnp')
			x_1.append('PON')
			x_1.append('KRG_common')
			x_1.append('KRG_rare')
			x_1.append('ExAC_common')
			x_1.append('Germline_filter')
			x_1.append('WhileList')

			header = '\t'.join(x_1)
			wfile.write('%s\n' % header)
			wfile2.write('%s\n' % header)

		else:
			x_1 = x.replace('\n','').split('\t')
			chr = x_1[4]
			pos = x_1[5]
			ref = x_1[10]
			var = x_1[12]
			type = x_1[9]
			"""
			if type == 'DEL' or type == 'INS':
				pass
			else:

				dic_key = '%s_%s_%s_%s' % (chr, pos, ref, var)

				dbsnp_com = db_com(dic_key, dbsnp_dic)
				PON_com = db_com(dic_key, PON_dic)
				KRG_common_com = db_com(dic_key, KRG_common_dic)
				KRG_rare_com = db_com(dic_key, KRG_rare_dic)
				ExAC_com = exac(x_1[ExAC_AF_EAS_pos])

				x_1.append(dbsnp_com)
				x_1.append(PON_com)
				x_1.append(KRG_common_com)
				x_1.append(KRG_rare_com)
				x_1.append(ExAC_com)

				if dbsnp_com == 'Yes' or PON_com == 'Yes' or KRG_common_com == 'Yes' or KRG_rare_com == 'Yes' or ExAC_com == 'Yes':
					x_1.append('Yes')
					total_com = 'Yes'
				else:
					x_1.append('No')
					total_com = 'No'

				line = '\t'.join(x_1)
				wfile.write('%s\n' % line)
				
				if total_com == 'No' and exon_dic.has_key(x_1[8]):
					wfile2.write('%s\n' % line)
	wfile.close()
	"""



			dic_key = '%s_%s_%s_%s' % (chr, pos, ref, var)

			dbsnp_com = db_com(dic_key, dbsnp_dic)
			PON_com = db_com(dic_key, PON_dic)
			KRG_common_com = db_com(dic_key, KRG_common_dic)
			KRG_rare_com = db_com(dic_key, KRG_rare_dic)
			ExAC_com = exac(x_1[ExAC_AF_EAS_pos])

			x_1.append(dbsnp_com)
			x_1.append(PON_com)
			x_1.append(KRG_common_com)
			x_1.append(KRG_rare_com)
			x_1.append(ExAC_com)

			if dbsnp_com == 'Yes' or PON_com == 'Yes' or KRG_common_com == 'Yes' or KRG_rare_com == 'Yes' or ExAC_com == 'Yes':
				x_1.append('Yes')
				total_com = 'Yes'
			else:
				x_1.append('No')
				total_com = 'No'
			x_1.append('No')

			line = '\t'.join(x_1)
			wfile.write('%s\n' % line)
			
			if total_com == 'No' and exon_dic.has_key(x_1[8]):
				wfile2.write('%s\n' % line)
	wfile.close()

if __name__=="__main__":__main__()
