#! /usr/bin/python
import optparse, os, sys, subprocess
from subprocess import *
from Bio.Seq import Seq
from utils import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input MAF file for filtering')
parser.add_option('-o', dest='output', help='output')
parser.add_option('-r', dest='hrd', help='HRD gene list')
parser.add_option('-m', dest='mmr', help='MMR gene list')
#parser.add_option('-D', dest='drug', help='drug file')
parser.add_option('-d', dest='diagnosis', help='diagnosis')
parser.add_option('-c', dest='ini',help='ini file')

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
	#hrd_genes = {'BRCA1':'1', 'BRCA2':'1', 'ATM':'1', 'BARD1':'1', 'BRIP1':'1', 'CHEK2':'1', 'NBN':'1', 'PALB2':'1', 'RAD51C':'1', 'MRE11A':'1', 'RAD51D':'1', 'ATR':'1', 'FAM175A':'1', 'XRCC2':'1'}
	hrd_genes = {}
	for x in open(options.hrd).xreadlines():
		x_1 = x.strip()
		hrd_genes[x_1] = '1'

	if hrd_genes.has_key(gene):
		val = 'Yes'
	else:
		val = 'No'
	return val


def mmr(gene):
	#mmr_genes = {'RAD51D':'1', 'RAD51C':'1', 'NBN':'1', 'XRCC2':'1', 'GEN1':'1', 'FAM175A':'1', 'MRE11A':'1', 'PMS2':'1', 'BRIP1':'1', 'PALB2':'1', 'BARD1':'1', 'MSH6':'1', 'MSH2':'1', 'BAP1':'1', 'BRCA1':'1', 'BRCA2':'1', 'MLH1':'1', 'ATR':'1', 'CHEK2':'1', 'ATM':'1'}
	mmr_genes = {}

	for y in open(options.mmr).xreadlines():
		y_1 = y.strip()
		mmr_genes[y_1] = '1'

	if mmr_genes.has_key(gene):
		val = 'Yes'
	else:
		val = 'No'
	return val

#def matching_drug(drug_db, gene, alteration, diagnosis):
def matching_drug(drug_db, negative, gene, alteration, diagnosis):
	dic = {}
	negative_genes = {}
	#db_path = '/data/D161740/Reference/drug/OncoKB_drug_20161021.txt'
	db_path = drug_db

	#for x in open(db_path):
	for x in open(db_path).xreadlines():
		x_1 = x.strip().split('\t')
		key = '%s_%s' % (x_1[0].upper(), x_1[1].upper())
		if not dic.has_key(key):
			dic[key] = {}
			if dic[key].has_key(x_1[2].upper()):
				pass
			else:
				dic[key][x_1[2].upper()] = x_1[4].upper()
		else:
			if dic[key].has_key(x_1[2].upper()):
				pass
			else:
				dic[key][x_1[2].upper()] = x_1[4].upper()


	diag = diagnosis.upper()

	for y in open(negative).xreadlines():
		y_1 = y.strip().split('\t')
		key2 = '%s_%s' % (y_1[0].upper(), y_1[1].upper())
		negative_genes[key2] = y_1

	match_key = '%s_%s' % (gene.upper(), alteration.upper())
	match_key2 = '%s_%s' % (gene.upper(), diag)

	if negative_genes.has_key(match_key2):
		drug = 'No response to TKI'
	else:
		if dic.has_key(match_key):
			if dic[match_key].has_key(diag):
				drug = dic[match_key][diag]

			else:
				drug = '-'
		else:
			drug = '-'

	return drug




def __main__():
	check_option(options.input, "You need to use '-i' option with input MAF file\n")

	input_maf = options.input

	params = read_config(options.ini)

	#wname = input_maf.replace('.maf', '.selected.maf') 
	wname = options.output
	wfile = open(wname,'w') 
	

	for x in open(input_maf).xreadlines(): 
		if x[0] == '#':
			pass

		elif x[0:11] == 'Hugo_Symbol':
			x_1 = x.replace('\n','').split('\t')

			len_con = len(x_1)


			#hugo_symbol = x_1.index('Gene')
			hugo_symbol = x_1.index('Hugo_Symbol')
			transcript = x_1.index('Transcript_ID')
			Variant_Type = x_1.index('Variant_Type')
			Chromosome = x_1.index('Chromosome')
			Start_Position = x_1.index('Start_Position')
			End_Position = x_1.index('End_Position')
			Reference_Allele = x_1.index('Reference_Allele')
			Tumor_Seq_Allele2 = x_1.index('Tumor_Seq_Allele2')

			Variant_Clasification = x_1.index('Variant_Classification')
			Variant_Type = x_1.index('Variant_Type')
			dbSNP = x_1.index('dbSNP_RS')
			mutation_status = x_1.index('Mutation_Status')

			Tumor_sample = x_1.index('Tumor_Sample_Barcode')
			Normal_sample = x_1.index('Matched_Norm_Sample_Barcode')

			Alteration_matched_with_db = x_1.index('Alteration_matched_with_db')

			exon = x_1.index('EXON')
			intron = x_1.index('INTRON')

			t_pos_total = x_1.index('t_depth')
			t_pos_ref = x_1.index('t_ref_count')
			t_pos_var = x_1.index('t_alt_count')
			n_pos_total = x_1.index('n_depth')
			n_pos_ref = x_1.index('n_ref_count')
			n_pos_var = x_1.index('n_alt_count')
			try:
				strand = x_1.index('STRAND')
			except:
				strand = x_1.index('STRAND_VEP')
			HGVSc = x_1.index('HGVSc')
			HGVSp = x_1.index('HGVSp_Short')
			#exon_number = x_1.index('Exon_Number')

			SIFT = x_1.index('SIFT')
			PolyPhen = x_1.index('PolyPhen')
			domain = x_1.index('DOMAINS')

			Variant_allele_ratio = x_1.index('Variant_allele_ratio')
			#Significant_call = x_1.index('Significant_call(VAF>=0.05)')
			#Tier = x_1.index('Tier')
			#Gene_codon = x_1.index('Gene_Codon_ProteinChange')
			#Pfam_domain = x_1.index('Pfam domain')
			white_list = x_1.index('WhiteList')
			#Sequence = x_1.index('Up-/Down-stream sequences')


			line = []
			#line.append(x_1[hugo_symbol])
			line.append('Gene')
			line.append('Transcript')
			line.append('Mutation_status')
			#line.append(x_1[Chromosome])
			line.append('Chr')
			line.append(x_1[Start_Position])
			line.append(x_1[End_Position])
			line.append(x_1[Reference_Allele])
			line.append(x_1[Tumor_Seq_Allele2])
			#line.append(x_1[Variant_Clasification])
			line.append('VarClass')
			#line.append(x_1[Variant_Type])
			line.append('Type')
			line.append(x_1[HGVSc])
			#line.append(x_1[HGVSp])
			line.append('Alteration')
			line.append('Alteration_matched_with_db')
			#line.append(x_1[exon_number])
			#line.append(x_1[Gene_codon])
			#line.append(x_1[Tier])
			line.append(x_1[dbSNP])
			line.append(x_1[t_pos_total])
			line.append(x_1[t_pos_ref])
			line.append(x_1[t_pos_var])
			line.append(x_1[n_pos_total])
			line.append(x_1[n_pos_ref])
			line.append(x_1[n_pos_var])
			#line.append(x_1[Variant_allele_ratio])
			line.append('Allele frequency')
			#line.append(x_1[Significant_call])
			line.append(x_1[strand])
			line.append(x_1[exon])
			line.append(x_1[intron])
			line.append(x_1[SIFT])
			line.append(x_1[PolyPhen])
			line.append(x_1[domain])
			#line.append(x_1[Pfam_domain])
			#line.append(x_1[Sequence])
			line.append('HRD')
			line.append('MMR')
			line.append('Zygosity')
			line.append('TranscriptRank')
			line.append('Diagnosis')
			line.append('Drug FDA Approved (same tumor type)')
			line.append('Drugable')
			line.append('WhiteList')


			cmd = '\t'.join(line)
			wfile.write('%s\n' % cmd)

		else:


			x_1 = x.replace('\n','').split('\t')

			line = []
			line.append(x_1[hugo_symbol])
			line.append(x_1[transcript])
			line.append(x_1[mutation_status])
			line.append(x_1[Chromosome])
			line.append(x_1[Start_Position])
			line.append(x_1[End_Position])
			line.append(x_1[Reference_Allele])
			line.append(x_1[Tumor_Seq_Allele2])
			line.append(x_1[Variant_Clasification])
			#line.append(x_1[Variant_Type])
			if x_1[Variant_Type] == 'SNP':
				line.append('SNV')
			else:
				line.append(x_1[Variant_Type])
			line.append(x_1[HGVSc])
			line.append(x_1[HGVSp].replace('p.',''))
			line.append(x_1[Alteration_matched_with_db])
			#line.append(x_1[exon_number])
			#line.append(x_1[Gene_codon])
			#line.append(x_1[Tier])
			line.append(x_1[dbSNP])
			line.append(x_1[t_pos_total])
			line.append(x_1[t_pos_ref])
			line.append(x_1[t_pos_var])
			line.append(x_1[n_pos_total])
			line.append(x_1[n_pos_ref])
			line.append(x_1[n_pos_var])
			line.append(x_1[Variant_allele_ratio])
			#line.append(x_1[Significant_call])
			line.append(x_1[strand])
			line.append(x_1[exon])
			line.append(x_1[intron])
			line.append(x_1[SIFT])
			line.append(x_1[PolyPhen])
			line.append(x_1[domain])
			#line.append(x_1[Pfam_domain])
			#line.append(x_1[Sequence])
			line.append(hrd(x_1[hugo_symbol]))
			line.append(mmr(x_1[hugo_symbol]))
			line.append('Heterozygosity')
			line.append('1')
			line.append(options.diagnosis.replace('"',''))
			#line.append(matching_drug(params['actionable_variants'], params['negative_gene_list'], x_1[hugo_symbol], x_1[HGVSp].replace('p.',''), options.diagnosis.replace('"','')))
			line.append(matching_drug(params['actionable_variants'], params['negative_gene_list'], x_1[hugo_symbol], x_1[Alteration_matched_with_db], options.diagnosis.replace('"','')))
			#if matching_drug(params['actionable_variants'], params['negative_gene_list'], x_1[hugo_symbol], x_1[HGVSp].replace('p.',''), options.diagnosis.replace('"','')) == '-':
			if matching_drug(params['actionable_variants'], params['negative_gene_list'], x_1[hugo_symbol], x_1[Alteration_matched_with_db], options.diagnosis.replace('"','')) == '-':
				line.append('No')
			else:
				if x_1[mutation_status] == 'Somatic':
					line.append('Yes')
				else:
					line.append('No')
			line.append(x_1[white_list])


			cmd = '\t'.join(line)
			wfile.write('%s\n' % cmd)



	wfile.close()
if __name__=="__main__":__main__()
