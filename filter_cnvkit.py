#! /usr/bin/python
import optparse, os, sys, subprocess, re
from subprocess import *
from Bio.Seq import Seq
import shutil

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input MAF file for filtering')

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

def hrd(gene):
	hrd_genes = {'BRCA1':'1', 'BRCA2':'1', 'ATM':'1', 'BARD1':'1', 'BRIP1':'1', 'CHEK2':'1', 'NBN':'1', 'PALB2':'1', 'RAD51C':'1', 'MRE11A':'1', 'RAD51D':'1', 'ATR':'1', 'FAM175A':'1', 'XRCC2':'1'}
	if hrd_genes.has_key(gene):
		val = 'Yes'
	else:
		val = 'No'
	return val


def mmr(gene):
	mmr_genes = {'RAD51D':'1', 'RAD51C':'1', 'NBN':'1', 'XRCC2':'1', 'GEN1':'1', 'FAM175A':'1', 'MRE11A':'1', 'PMS2':'1', 'BRIP1':'1', 'PALB2':'1', 'BARD1':'1', 'MSH6':'1', 'MSH2':'1', 'BAP1':'1', 'BRCA1':'1', 'BRCA2':'1', 'MLH1':'1', 'ATR':'1', 'CHEK2':'1', 'ATM':'1'}
	if mmr_genes.has_key(gene):
		val = 'Yes'
	else:
		val = 'No'
	return val

def sigflag(gene, cn):
	drug = {'ERBB2':'Amplification', 'CDK4':'Amplification', 'MET':'Amplification', 'FGFR1':'Amplification', 'FGFR2':'Amplification'}
	#if cn >= 4:
	if cn >= 5:
		event = "Amplifictiona"
	elif cn == 0:
		event = "Deletion"
	
	if drug.has_key(gene) and drug[gene] == event:
		val = 'Yes'
	else:
		val = 'No'
	return val

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

"""
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
"""




raw_cns = options.input.replace('.call.cns', '.cns')

wname = options.input.replace('.cns', '.fixed.cns')
wname2 = options.input.replace('.cns', '.fixed.splited.cns')
wfile = open(wname, 'w')
wfile2 = open(wname2, 'w')

black = {'TRG':'', 'TRB':'', 'TRA':'', 'IGK':'', 'IGH':'', 'IGL':''}

drug = {'ERBB2':'Amplification', 'CDK4':'Amplification', 'MET':'Amplification', 'FGFR1':'Amplification', 'FGFR2':'Amplification'}

cns_dic = {}

for k in open(raw_cns).xreadlines():
	k_1 = k.strip().split('\t')
	cns_key = '%s_%s_%s' % (k_1[0], k_1[1], k_1[2])
	cns_dic[cns_key] = k_1[4]


for x in open(options.input).xreadlines():
	x_1 = x.strip().split('\t')
	key2 = '%s_%s_%s' % (x_1[0], x_1[1], x_1[2])

	if x_1[0] == 'chromosome':
		wfile.write(x)
		header = 'Type\tChr\tStart\tEnd\tGene\tlog2\tCN\tAlteration\tSigFlag\tMMR\tHRD\tCnvType\tKnownFlag\n'
		wfile2.write(header)

	else:
		genes = x_1[3].split(',')
		#cn = int(x_1[7])
		cn = int(x_1[5])

		if x_1[3] != '-':
			if cn == 0 or cn >= 4:
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
							#con = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('CNV', x_1[0], x_1[1], x_1[2], y, x_1[4], x_1[5], copy_num(cn), sigflag(y, cn), mmr(y), hrd(y), 'CNV', sigflag(y, cn))
							con = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('CNV', x_1[0], x_1[1], x_1[2], y, cns_dic[key2], x_1[5], copy_num(cn), sigflag(y, cn), mmr(y), hrd(y), 'CNV', sigflag(y, cn))
							wfile2.write(con)
			else:
				if black.has_key(x_1[3]):
					pass
				else:
					if copy_num(cn) != 'Neutral':
						#con = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('CNV', x_1[0], x_1[1], x_1[2], x_1[3], x_1[4], x_1[5], copy_num(cn), sigflag(x_1[3], cn), mmr(x_1[3]), hrd(x_1[3]), 'CNV', sigflag(x_1[3], cn))
						con = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('CNV', x_1[0], x_1[1], x_1[2], x_1[3], cns_dic[key2], x_1[5], copy_num(cn), sigflag(x_1[3], cn), mmr(x_1[3]), hrd(x_1[3]), 'CNV', sigflag(x_1[3], cn))
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


