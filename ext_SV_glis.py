#! /usr/bin/python
import optparse, os, sys, subprocess, re, MySQLdb, glob, shutil

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-d', dest='dir', help='analyzed dir e.g) 131016_M02071_0004_000000000-A50PV')
parser.add_option('-o', dest='output', help='output')
parser.add_option('-r', dest='hrd', help='HRD gene list')
parser.add_option('-m', dest='mmr', help='MMR gene list')
parser.add_option('-t', dest='target', help='SV target gene')
parser.add_option('-p', dest='pon', help='breakmer PON file')


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
		mmr_genes[y_1] ='1'

	if mmr_genes.has_key(gene):
		val = 'Yes'
	else:
		val = 'No'
	return val

def target(gene):
	genes = gene.split('-')

	#target = {'NTRK1':'', 'ALK':'', 'ROS1':'', 'EGFR':'', 'ABL1':'', 'RET':'', 'TMPRSS2':'', 'EWSR1':''}
	target = {}

	for z in open(options.target).xreadlines():
		z_1 = z.strip()
		target[z_1] = '1'

	val = ''

	for x in genes:
		if target.has_key(x):
			val = 'Yes'
	if val == '':
		val = '-'
	return val

def __main__():
	check_option(options.dir, "You need to use '-d' option with basecall dir\n")

	target_dir = '%s/output' % options.dir

	#files = glob.glob('%s/*_svs.out' % target_dir)
	files = glob.glob('%s/*_rearrangement_svs.out' % target_dir)

	if len(files) == 0:
		wname = options.output
		wfile = open(wname, 'w')
		header = 'Gene\ttarget_breakpoints\talign_cigar\tmismatches\tstrands\trep_overlap_segment_len\tsv_type\tsplit_read_count\tnkmers\tdisc_read_count\tbreakpoint_coverages\tcontig_id\tcontig_seq\tType\tAlteration\tAlteration_matched_with_db\tTransType\tRearrangement_target\n'
		wfile.write(header)
		wfile.close()
		sys.exit()

	else:

		#pon_file = '/data/D161740/Reference/PON/SV/breakmer_PON_20151029.txt'
		#pon_file = '/data/D161740/Reference/PON/SV/breakmer_PON_20160926.txt'
		pon_file = options.pon

		dic = {}

		for z in open(pon_file).xreadlines():
			z_1 = z.strip().split('\t')
			key = '%s\t%s' % (z_1[0], z_1[1])
			dic[key] = 0


		header = []
		con = []
		
		a = 0

		for x in files:
			for no, y in enumerate(open(x).readlines()):
				if a == 0 and no == 0:
					if len(header) == 0:
						header = 'Gene\ttarget_breakpoints\talign_cigar\tmismatches\tstrands\trep_overlap_segment_len\tsv_type\tsplit_read_count\tnkmers\tdisc_read_count\tbreakpoint_coverages\tcontig_id\tcontig_seq\tType\tAlteration\tAlteration_matched_with_db\tTransType\tRearrangement_target\n'
						con.append(header)
						a += 1

				if no != 0:
					a += 1
					y_1 = y.strip().split('\t')
					m_key = '%s\t%s' % (y_1[0], y_1[1])
					if dic.has_key(m_key):
						pass
					else:
						if y_1[6] == 'rearrangement':
							y_1[0] = y_1[0].replace(',','-')
							y2 = '\t'.join(y_1)
							type = 'Rearrangement'
							con.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (y2, type, 'Fusion', 'Fusion', 'Intra', target(y_1[0])))


		con_len = len(con)
		print con_len
		contents = '\n'.join(con)


		if con_len != 1:
			wname = options.output
			wfile = open(wname, 'w')
			wfile.write(contents)
			wfile.close()


if __name__=="__main__":__main__()
