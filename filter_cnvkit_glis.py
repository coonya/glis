#! /usr/bin/python
import optparse, os, sys, subprocess, re
from subprocess import *
from Bio.Seq import Seq
import shutil
from utils import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input MAF file for filtering')
parser.add_option('-r', dest='hrd',help='hrd list')
#parser.add_option('-D', dest='drug',help='drug list')
parser.add_option('-d', dest='diagnosis',help='diagnosis')
parser.add_option('-c', dest='ini',help='ini file')

(options, args) = parser.parse_args()

def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())


def read_list(db):
	dic = {}
	print db

	for x in open(db).xreadlines():
		x_1 = x.strip()
		dic[x_1] = 1
	return dic

def db_com(dic, gene, cn):
	key1 = '%s_%s' % (gene, cn)
	if dic.has_key(key1):
		val = 'Yes'
	else:
		val = 'No'
	return val

def db_com2(dic, gene):
	if dic.has_key(gene):
		val = 'Yes'
	else:
		val = 'No'
	return val


def drug_list(db, diagnosis):
	dic = {}
	level = {'1':'1', '2A':'1', '2B':'1'}

	for x in open(db).xreadlines():
		x_1 = x.strip().split('\t')
		#if x_1[1] == 'Amplification' and dic.has_key(x_1[3]) and diagnosis == x_1[2]:
		if x_1[1] == 'Amplification' and diagnosis == x_1[2]:
			key2 = '%s_%s' % (x_1[0], x_1[1])
			dic[key2] = x
	return dic


def copy_num(cn):
	#if cn >= 4:
	if cn >= 5:
		val = 'Amplification'
	elif cn == 0:
		val = 'Loss'
	else:
		val = 'Neutral'
	return val


check_option(options.input, "you need to use '-i' options.")

params = read_config(options.ini)

drug_db = params['actionable_variants']

raw_cns = options.input.replace('.call.cns', '.cns')

wname = options.input.replace('.cns', '.fixed.cns')
wname2 = options.input.replace('.cns', '.fixed.splited.cns')
wfile = open(wname, 'w')
wfile2 = open(wname2, 'w')

black = {'TRG':'', 'TRB':'', 'TRA':'', 'IGK':'', 'IGH':'', 'IGL':''}

cns_dic = {}
print options.hrd
print options.input
print options.diagnosis

hrd_dic = read_list(options.hrd)
drug_dic = drug_list(drug_db, options.diagnosis)

for k in open(raw_cns).xreadlines():
	k_1 = k.strip().split('\t')
	cns_key = '%s_%s_%s' % (k_1[0], k_1[1], k_1[2])
	cns_dic[cns_key] = k_1[4]


print cns_dic
for x in open(options.input).xreadlines():
	x_1 = x.strip().split('\t')
	key2 = '%s_%s_%s' % (x_1[0], x_1[1], x_1[2])

	if x_1[0] == 'chromosome':
		wfile.write(x)
		header = 'Type\tChr\tStart\tEnd\tGene\tlog2\tCN\tAlteration\tAlteration_matched_with_db\tSigFlag\tHRD\tCnvType\tKnownFlag\n'
		wfile2.write(header)

	else:
		genes = x_1[3].split(',')
		#cn = int(x_1[7])
		cn = int(x_1[5])

		if x_1[3] != '-':
			#if cn == 0 or cn >= 4:
			if cn == 0 or cn >= 5:
				x_1[4] = cns_dic[key2]
				line = '\t'.join(x_1)
				wfile.write('%s\n' % line)
			else:
				pass

			if len(genes) != 1:
				for y in genes:
					if black.has_key(y):
						pass
					else:
						if copy_num(cn) != 'Neutral':
							x_1[3] = y
							if x_1[0] != 'X':
							#con = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('CNV', x_1[0], x_1[1], x_1[2], y, x_1[4], x_1[5], copy_num(cn), sigflag(y, cn), mmr(y), hrd(y), 'CNV', sigflag(y, cn))
								con = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('CNV', x_1[0], x_1[1], x_1[2], y, cns_dic[key2], x_1[5], copy_num(cn), copy_num(cn), db_com(drug_dic, y, copy_num(cn)), db_com2(hrd_dic, y), 'CNV', db_com(drug_dic, y, copy_num(cn)))
								wfile2.write(con)
			else:
				if black.has_key(x_1[3]):
					pass
				else:
					if copy_num(cn) != 'Neutral':
						if x_1[0] != 'X':
						#con = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('CNV', x_1[0], x_1[1], x_1[2], x_1[3], x_1[4], x_1[5], copy_num(cn), sigflag(x_1[3], cn), mmr(x_1[3]), hrd(x_1[3]), 'CNV', sigflag(x_1[3], cn))
							con = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('CNV', x_1[0], x_1[1], x_1[2], x_1[3], cns_dic[key2], x_1[5], copy_num(cn), copy_num(cn),db_com(drug_dic, x_1[3], copy_num(cn)), db_com2(hrd_dic, x_1[3]), 'CNV', db_com(drug_dic, x_1[3], copy_num(cn)))
							wfile2.write(con)


wfile.close()
wfile2.close()




### change file name
if '.fastq.gz' in options.input:
	template = '.fastq.gz.initialAlign.merged.dedup.realign'

	cnr = options.input.replace('.call.cns', '.cnr')
	filtered_cnr = options.input.replace('.call.cns', '.filtered.cnr')
	cns = options.input.replace('.call.cns', '.cns')
	call_cns = options.input
	fixed_cns = options.input.replace('.call.cns', '.call.fixed.cns')
	splited_cns = options.input.replace('.call.cns', '.call.fixed.splited.cns')

	shutil.copy2(cnr, cnr.replace(template, ''))
	shutil.copy2(filtered_cnr, filtered_cnr.replace(template, ''))
	shutil.copy2(cns, cns.replace(template, ''))
	shutil.copy2(call_cns, call_cns.replace(template, ''))
	shutil.copy2(fixed_cns, fixed_cns.replace(template,''))
	shutil.copy2(splited_cns, splited_cns.replace(template,''))

