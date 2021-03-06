#! /usr/bin/python
import optparse, os, sys, subprocess
from subprocess import *
from Bio.Seq import Seq
from utils import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input MAF file for filtering')
parser.add_option('-c', dest='ini',help='ini file')
#parser.add_option('-c', dest='caller',help='somatic mutation caller e.g.) mutect, somaticindelocator')

(options, args) = parser.parse_args()

def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())



def check_region(line):
	regions = {"3'Flank": "",\
			"3'UTR": "",\
			"Frame_Shift_Del": "",\
			"Frame_Shift_Ins": "",\
			"IGR": "",\
			"In_Frame_Del": "",\
			"In_Frame_Ins": "",\
			"Intron": "",\
			"RNA": "",\
			"Missense_Mutation": "",\
			"Nonsense_Mutation": "",\
			"Silent": "",\
			"Splice_Site": ""}
	exonic = {"Frame_Shift_Del": "",\
			"Frame_Shift_Ins": "",\
			"In_Frame_Del": "",\
			"In_Frame_Ins": "",\
			"Missense_Mutation": "",\
			"Nonsense_Mutation": "",\
			"Nonstop_Mutation": "",\
			"Splice_Site": ""}

	region = line.strip().split('\t')[8]
	if exonic.has_key(region):
		return 1
	else:
		return 0

def reverse_complement(seq):
	seq = Seq(seq)
	return str(seq.reverse_complement())


def significance(t_count_ref, t_count_var, var_ratio):
	if (int(t_count_ref) + int(t_count_var)) >= 30 and float(var_ratio) >= 0.05:
		#if float(var_ratio) >= 0.05:
		val = 'Yes'
	else:
		val = 'No'
	
	return val


def determine_type(anno_dic, action_dic, pos_dic, line):
	x_1 = line.replace('\n', '').split('\t')
	print x_1
	

	var = x_1[pos_dic['protein']].split('.')[1]
	key_12 = '%s_%s' % (x_1[0], var)
	if action_dic.has_key(key_12):
		res_f = var
	elif check_region(line) == 1:
		if x_1[pos_dic['var_class']] == 'Nonsense_Mutation':
			res = 'Truncating Mutations'
			key_13 = '%s_%s' % (x_1[0], res)
			if action_dic.has_key(key_13):
				res_f = res
			else:
				res_f = var
		#elif x_1[pos_dic['Variant_Type']] == 'DEL' or x_1[pos_dic['Variant_Type']] == 'INS':
		else:
			if x_1[pos_dic['Variant_Type']] == 'DEL':
				val = 'deletion'
			elif x_1[pos_dic['Variant_Type']] == 'INS':
				val = 'insertion'
			else:
				val = 'mutations'

			val2 = 'Exon %s %s' % ( x_1[pos_dic['exon_number']].split('/')[0], val)
			key_13 = '%s_%s' % (x_1[0], val2)
			if action_dic.has_key(key_13):
				res_f = val2
			else:
				res_f = var
	else:
		if anno_dic.has_key(key_12):
			can_anno = anno_dic[key_12].split('\t')
			if can_anno[2] == 'Oncogenic':
				val = 'Oncogenic Mutations'
				if action_dic.has_key(val):
					res_f = val
				else:
					res_f = var
		else:
			res_f = var
	print res_f
	return res_f
			
	



def __main__():
	check_option(options.input, "You need to use '-i' option with input MAF file\n")

	input_maf = options.input

	params = read_config(options.ini)
	
	anno_dic = {}
	
	for anno in open(params['annotated_variants']).xreadlines():
		anno_1 = anno.strip().split('\t')
		key11 = '%s_%s' % (anno_1[0], anno_1[1])
		anno_dic[key11] = anno
	
	action_dic = {}

	for k in open(params['actionable_variants']).xreadlines():
		k_1 = k.strip().split('\t')
		key22 = '%s_%s' % (k_1[0], k_1[1])
		action_dic[key22] = k
	
	wname = input_maf.replace('.maf', '.exon.maf') 
	wname2 = input_maf.replace('.maf', '.all.maf') 
	wfile = open(wname,'w') 
	wfile2 = open(wname2,'w') 
	
	pos_dic = {}
	

        for x in open(input_maf).xreadlines(): 
		if x[0] == '#':
			pass

		elif x[0:11] == 'Hugo_Symbol':
			x_1 = x.replace('\n','').split('\t')

			len_con = len(x_1)

			ens_gene_pos = x_1.index('Gene')
			x_1[ens_gene_pos] = 'ENS_gene'

			hugo_symbol = x_1.index('Hugo_Symbol')
			Variant_Type = x_1.index('Variant_Type')
			chromosome = x_1.index('Chromosome')
			Start_Position = x_1.index('Start_Position')
			Reference_Allele = x_1.index('Reference_Allele')
			Tumor_Seq_Allele2 = x_1.index('Tumor_Seq_Allele2')

			center = x_1.index('Center')
			exon_number = x_1.index('Exon_Number')
			exon = x_1.index('EXON')
			intron = x_1.index('INTRON')
			var_class = x_1.index('Variant_Classification')

			t_pos_total = x_1.index('t_depth')
			t_pos_ref = x_1.index('t_ref_count')
			t_pos_var = x_1.index('t_alt_count')
			n_pos_total = x_1.index('n_depth')
			n_pos_ref = x_1.index('n_ref_count')
			n_pos_var = x_1.index('n_alt_count')
			strand = x_1.index('STRAND_VEP')

			mutation_status = x_1.index('Mutation_Status')
			#sequence_source = x_1.index('Sequence_Source')
			sequencer = x_1.index('Sequencer')
			codon = x_1.index('HGVSc')
			protein = x_1.index('HGVSp_Short')
			cdna_pos = x_1.index('CDS_position')
			domain = x_1.index('DOMAINS')

			x_1.append('Variant_allele_ratio')
			x_1.append('Significant_call(VAF>=0.05)')
			x_1.append('Gene_Codon_ProteinChange')
			x_1.append('TranscriptRank')
			x_1.append('Alteration_matched_with_db')
			cmd = '\t'.join(x_1)
			wfile.write('%s\n' % cmd)
			wfile2.write('%s\n' % cmd)
			pos_dic['protein'] = protein
			pos_dic['var_class'] = var_class
			pos_dic['Variant_Type'] = Variant_Type
			pos_dic['exon_number'] = exon_number

		else:
			exonic = check_region(x)

			x_1 = x.replace('\n','').split('\t')

			if len(x_1) != len_con:
				for le in range(len(x_1), len_con):
					x_1.append('-')

			for no, k in enumerate(x_1):
				if k == '':
					x_1[no] = '-'
			gene = x_1[hugo_symbol]
			mutation_type = x_1[Variant_Type]
			chr = x_1[chromosome]

			if mutation_type == 'DEL':
				pos = str(int(x_1[Start_Position])-1)
			else:
				pos = x_1[Start_Position]

			ref = x_1[Reference_Allele]
			alt = x_1[Tumor_Seq_Allele2]

			if mutation_type == 'DEL':
				exonic_key = '%s_%s' % (x_1[chromosome], int(x_1[Start_Position])-1)
			else: 
				exonic_key = '%s_%s' % (x_1[chromosome], x_1[Start_Position])


			variant1 = '%s' % (exonic_key)
			variant2 = '%s_%s_%s' % (exonic_key, ref, alt)

			
			x_1[center] = 'AMC'
			x_1[exon_number] = x_1[exon_number].replace('/','//')
			try:
				x_1[exon] = x_1[exon].replace('/','//')
			except:
				pass
			try:
				x_1[intron] = x_1[intron].replace('/','//')
			except:
				pass


			#x_1[mutation_status] = 'Somatic'
			#x_1[sequence_source] = panel(options.panel)


			x_1[sequencer] = 'Illumina MiSEQ'
			
			if x_1[t_pos_var] != '' and x_1[t_pos_total] != '':
				if x_1[t_pos_ref] == '-':
					x_1[t_pos_ref] = str(int(x_1[t_pos_total]) - int(x_1[t_pos_var]))
			
			if x_1[t_pos_var] == '-' or x_1[t_pos_total] == '-':
				var_ratio = 0
			else:
				var_ratio = '%5.2f' % float(float(x_1[t_pos_var]) / float(x_1[t_pos_total]))

			x_1.append(str(var_ratio))
			x_1.append(significance(x_1[t_pos_ref], x_1[t_pos_var], var_ratio))

			## fix codon
			if mutation_type == 'DEL' or mutation_type == 'INS':
				new_codon_final = x_1[codon]
			else:
				if mutation_type == 'SNP':
					if x_1[codon] != '':
						new_codon_pos = x_1[cdna_pos].split('/')[0]

						if x_1[strand] == '1':
							new_codon_ref = x_1[Reference_Allele]
							new_codon_alt = x_1[Tumor_Seq_Allele2]
						elif x_1[strand] == '-1':
							new_codon_ref = reverse_complement(x_1[Reference_Allele])
							new_codon_alt = reverse_complement(x_1[Tumor_Seq_Allele2])
						else:
							new_codon_ref = x_1[Reference_Allele]
							new_codon_alt = x_1[Tumor_Seq_Allele2]

						new_codon_final = 'c.%s%s>%s' % (new_codon_pos, new_codon_ref, new_codon_alt)
						x_1[codon] = new_codon_final

					else:
						if x_1[strand] == '1':
							new_codon = list(str(x_1[codon]))
							new_codon[-3] = x_1[Reference_Allele]
							new_codon[-1] = x_1[Tumor_Seq_Allele2]
							new_codon_final = ''.join(new_codon)
						elif x_1[strand] == '-1':
							new_codon = list(str(x_1[codon]))
							new_codon[-3] = reverse_complement(x_1[Reference_Allele])
							new_codon[-1] = reverse_complement(x_1[Tumor_Seq_Allele2])
							new_codon_final = ''.join(new_codon)
				elif mutation_type =='DEL' or mutation_type == 'INS':
					new_codon_final = x_1[codon]

				else:
					new_codon = x_1[codon].split('del')
					new_codon_pos = new_codon[0]
					new_codon_ref = new_codon[1].split('ins')[0]
					new_codon_var = new_codon[1].split('ins')[1]

					if x_1[strand] == '1':
						new_codon_final = '%sdel%sins%s' % (new_codon_pos, x_1[Reference_Allele], x_1[Tumor_Seq_Allele2])
					elif x_1[strand] == '-1':
						new_codon_final = '%sdel%sins%s' % (new_codon_pos, reverse_complement(x_1[Reference_Allele]), reverse_complement(x_1[Tumor_Seq_Allele2]))


			## Gene_Codon_ProteinChange
			x_1.append('%s_%s_%s' % (gene, new_codon_final, x_1[protein]))	
			## transcript Rank
			x_1.append('1')

			#x_1.append(determine_type(anno_dic, action_dic, pos_dic, x))



			if exonic == 1:	
				x_1.append(determine_type(anno_dic, action_dic, pos_dic, x))
				cmd = '\t'.join(x_1)
				wfile.write('%s\n' % cmd)
			else:
				if x_1[protein] == '-':
					x_1.append('-')
				else:
					x_1.append(x_1[protein].split('.')[1])
				cmd = '\t'.join(x_1)
			wfile2.write('%s\n' % cmd)
	wfile.close()
	wfile2.close()
if __name__=="__main__":__main__()
