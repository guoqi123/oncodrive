import os
import sys
from optparse import OptionParser

from auxiliar import *

###################################################3
##
##   Parser 
###################################################3
def parse_options():
	parser = OptionParser()

	parser.add_option("-s", "--syn", dest="syn_mut_fp", type="string",
					  action="store",
					  help="The path to the Synonimous mutations file to construct the background model",
					  metavar="PATH")

	parser.add_option("-n", "--nonsyn", dest="non_syn_mut_fp", type="string",
					  action="store",
					  help="The path to the NON-Synonimous mutations file to be checked", metavar="PATH")

	parser.add_option("-p", "--pos", dest="pos_pos", type="int", default="-1",
					  action="store",
					  help="Position of the mut files in which the AA position is stated ('-1' by default)",
					  metavar="INT")

	parser.add_option("-m", "--muts", dest="min_gene_mutations", type="int", default="3",
					  action="store",
					  help="Minimum number of mutations of a gene to be included in the analysis ('3' by default)",
					  metavar="INT")

	parser.add_option("-o", "--out", dest="output_fp", type="string",
					  action="store",
					  help="States the output file path", metavar="PATH")

	options = parser.parse_args()[0]

	if not options.syn_mut_fp or not options.non_syn_mut_fp or not os.path.exists(
			options.syn_mut_fp) or not os.path.exists(options.non_syn_mut_fp):
		sys.exit(
			"\nERROR!\nPlease specify a valid path for both Syn and Nonsyn mutations files  by using the corresponding"
			" launcher arguments. Use --help for further details.\n")
	if not options.output_fp:
		sys.exit("\nERROR!\nPlease specify a path for output file by using the corresponding"
				 " launcher arguments. Use --help for further details.\n")
	return options


###################################################3
##
##   Output 
###################################################3

def create_output_file(options, non_syn_accum_mut_pos_dict, non_syn_cluster_coordinates_dict, \
					   syn_gene_cluster_scores_dict, non_syn_gene_cluster_scores_dict, \
					   non_syn_cluster_muts_dict, non_syn_gene_cluster_scores_external_z_dict):
	output_fp = options.output_fp
	min_gene_mutations = options.min_gene_mutations

	cgc_dict = f2dict(cancer_census_fp, 0, 1)

	#Output file
	f_out_fp = output_fp + '_scoreCorrection' + str(CLUSTER_SCORE_CORRECTION) + '_minMuts' + str(
		min_gene_mutations) + 'OncoCLUST.txt'
	f_out = open(f_out_fp, 'w')
	f_out.write('\t'.join(
		['Gene', 'CGC', 'Gene_len', 'Gene_Muts', 'n_clusters', 'Muts_in_clusters', 'Gene_cluster_score', 'Z_value',
		 'P_value', 'Q_value']))

	m_out = []
	m_not_included_out = []
	for gene in non_syn_gene_cluster_scores_dict:
		cgc = cgc_dict[gene] if gene in cgc_dict else ''
		gene_len = int(cds_dict[gene]) / 3
		gene_muts = sum([non_syn_accum_mut_pos_dict[gene][pos] for pos in non_syn_accum_mut_pos_dict[gene].keys()])
		n_clusters = len(non_syn_cluster_coordinates_dict[gene].keys())

		if n_clusters > 0:
			muts_clusters = sum(non_syn_cluster_muts_dict[gene].values())
			gene_additive_cluster_score = non_syn_gene_cluster_scores_dict[gene]
			z_ext = non_syn_gene_cluster_scores_external_z_dict[
				gene] if gene in non_syn_gene_cluster_scores_external_z_dict else 'NA'
			m_out.append(
				[gene, cgc, gene_len, gene_muts, n_clusters, muts_clusters, gene_additive_cluster_score, z_ext])
		else:
			m_not_included_out.append([gene, cgc, gene_len, gene_muts, 0, 0, 0, 'NA', 'NA', 'NA'])

	m_out = order_matrix(m_out, -1)
	m_out = add_fdr(m_out, -1)
	m_out += m_not_included_out

	for l_out in m_out:
		f_out.write('\n' + '\t'.join(str(e) for e in l_out))
	f_out.close()
	print '\n\n====File created {0}=====\n'.format(f_out_fp)
	return 


###################################################3
##
##   Called Functions 
###################################################3

'''
First, it calculate the meaningful positions which are those with the minimum number of positions that 
are expected by less than 1% of probability according to the binomial model
   -->	input : dict = {pos: num_mutations}  and gene_len, gene_muts
   -->  output: dict = {cluster_id: [positions]}
'''
def get_gene_clusters(gene, gene_len, gene_muts, gene_accum_mut_pos_dict):
	gene_cluster_dict = {}

	#The 'meaninful positions are those that have less than 1% of being found by chance. This is a binomial
	# cumulated probability model in which p(X >=n) = 1 - p(X <=(n-1)), where the parameters are: prob_of_success = 1/len_gene
	# and num_of_success = n, and num_of_trials = num_of_mutations

	#call func to return the minimum number of mutations per position to be included in the cluster
	# (note that 'NOISE_BINOM_P_THRESHOLD' is a global variable)
	minimum_mut_per_position_threshold = get_binomial_minimum_mut_per_position_threshold(gene_len, gene_muts, NOISE_BINOM_P_THRESHOLD, MIN_MUTS_PER_MEANINGF_POS)

	#given the minimum # of muts per position to be included in a cluster, retrieve the positions
	significant_gene_positions_l = get_significant_gene_positions(gene_accum_mut_pos_dict, minimum_mut_per_position_threshold)

	#if no meaningful positions, no cluster can be obtained  
	if not len(significant_gene_positions_l) > 0: 
		return {}

	#given the significant positions, group them
	# (applying the constraint of NOISE_BINOM_P_THRESHOLD, which is a global variable)
	# gene_dict = {cluster_id: [meaningful_pos_init, meaningful_pos_end]}
	gene_meaningful_cluster_dict = group_meaningful_positions(gene, significant_gene_positions_l)

	#given the meaningful positions, check if there is some non meaningful position around
	gene_extended_meaningful_cluster_dict = extend_meaningful_cluster(gene, gene_meaningful_cluster_dict, gene_accum_mut_pos_dict)

	#print 'Gene {0} and cluster dict is {1}'.format(gene, gene_meaningful_cluster_dict)
	#print 'Gene {0} and cluster extended dict is {1}'.format(gene, gene_extended_meaningful_cluster_dict)

	return gene_extended_meaningful_cluster_dict

'''
returns a dict with the number (absolute) of mutations
in each position of each gene: {gene: {pos:accumulated_mutations}}
'''
def pos_f2accumulated_pos_dict(mut_fp, pos_pos):
	accum_pos_dict = {}
	pos_f = open(mut_fp, 'r')
	pos_f.readline()
	for l in pos_f:
		l_s = l.rstrip().split('\t')
		try:
			gene, pos = l_s[0], int(l_s[pos_pos])
			if gene not in accum_pos_dict:
				accum_pos_dict[gene] = {pos: 1}
			elif pos not in accum_pos_dict[gene]:
				accum_pos_dict[gene][pos] = 1
			else:
				accum_pos_dict[gene][pos] += 1
		except:
			#the position entry contains a non valid value (i.e. a non integer)
			continue
	pos_f.close()
	return accum_pos_dict


'''
(A) Define meaningful positions in each gene, i.e. those positions showing a non-syn mutation occurrence higher than expected by chance.

(B) The meaningful positions are clustered into groups that respect a maximal distance between two meaningful positions of the cluster

Note that the input is a dict with the number of mutations per position: {gene: {pos:accumulated_mutations}}

And returns a duct with the cluster_coordinates per gene (and per cluster): {gene: {cluster_id: [init_position, end_position]}}
'''
def cluster_mutations(accum_mut_pos_dict, min_gene_mutations):
	genes_sorted = sorted(accum_mut_pos_dict.keys())

	cluster_coordinates_dict = {}
	#---> Retrieval of clusters per gene
	for gene in genes_sorted:
		#---> (1) get the number of mutations of the gene 
		gene_muts = sum([accum_mut_pos_dict[gene][pos] for pos in accum_mut_pos_dict[gene].keys()]) 

		# if there are few muts, as stated by the user (or default threshold),the gene is not included
		if gene_muts < min_gene_mutations:
			continue
		#I have not the info of the gene len, ,the gene is not included 
		if gene not in cds_dict:
			continue

		gene_len = int(cds_dict[gene]) / 3   #note that the cds is in bp, so i have to convert to AAs

		#---> (2) get the clusters of the gene, i.e. (a) remove positions with a observed mutation number below 99% of chance
		#         according to binomial cumulated probability; (b) the remaining 'meaningful' positions are grouped by
		#         using with the constraint that two positios of the same cluster can not be separated by more than 'n' AAs  
		#   dict = {gene: {cluster_id: [pos_init, pos_end]}}		
		cluster_coordinates_dict[gene] = get_gene_clusters(gene, gene_len, gene_muts, accum_mut_pos_dict[gene])

	return cluster_coordinates_dict

'''
 inputs:
  #dict = {gene: {cluster_id: [init_position, end_position]}}
  #dict = {gene:{pos: num_mutations}}

And returns: dict = {gene: {cluster_id: num_enclosed_mutations}}
  (note that within the enclosed mutations, we take into account those that have been considered as being 
  in peak positions)
'''
def get_enclosed_mutations(cluster_coordinates_dict, accum_mut_pos_dict):
	cluster_summary_dict = {}
	for gene in cluster_coordinates_dict.keys():
		tmp_dict = {}
		for cluster_id in cluster_coordinates_dict[gene]:
			n_mutations = 0
			for position in accum_mut_pos_dict[gene]:
				cluster_init, cluster_end = cluster_coordinates_dict[gene][cluster_id][0],\
								 cluster_coordinates_dict[gene][cluster_id][1]
				if position in range(cluster_init, cluster_end +1):
					n_mutations += accum_mut_pos_dict[gene][position]
			tmp_dict[cluster_id] = n_mutations
		cluster_summary_dict[gene] = tmp_dict		
	return cluster_summary_dict


'''
CAlculate the scores for each cluster:
- input: {gene:{cluster_id: [pos_init, pos_end]}}, {gene:{pos: muts}}
- output: {gene: {cluster_id: score}}

'''
def score_clusters(accum_mut_pos_dict, cluster_coordinates_dict, cluster_muts_dict):
	cluster_scores_dict = {}
	genes_sorted = sorted(cluster_coordinates_dict.keys())
	for gene in genes_sorted:
		#---> (1) get the number of mutations of the gene 
		gene_muts = sum([accum_mut_pos_dict[gene][pos] for pos in accum_mut_pos_dict[gene].keys()]) 
		cluster_scores_dict[gene] = get_cluster_score_dict(gene, gene_muts, \
								accum_mut_pos_dict[gene], cluster_coordinates_dict[gene])
	return cluster_scores_dict

'''
The clustering score of a gene is the sum of the scores of all the clusters of that gene
'''
def group_clusters_of_genes(cluster_scores_dict):	
	gene_cluster_scores_dict = {}
	for gene in cluster_scores_dict:
			gene_cluster_scores_dict[gene] = sum(cluster_scores_dict[gene].values())
	return gene_cluster_scores_dict

'''
Inputs:
	dict = {gene: {cluster_id: [pos_init, pos_end]}}
	dict = {gene: {pos: num_mutations}}
The score here is calculated as the sum of the percentage of the mutations in each cluster pos, where each is divided
by the distance from the pos with maximal mutations
'''
def get_cluster_score_dict(gene, gene_mutations, gene_accum_pos_dict, gene_cluster_coordinates_dict):
	gene_clustering_score_dict = {}
	for cluster_id in gene_cluster_coordinates_dict:
		cl_init, cl_end = gene_cluster_coordinates_dict[cluster_id][0], gene_cluster_coordinates_dict[cluster_id][1]

		#look for the pos of the cluster with more mutations
		max_muts_pos = cl_init
		for pos in range(cl_init, cl_end + 1):
			if pos in gene_accum_pos_dict and gene_accum_pos_dict[pos] > gene_accum_pos_dict[max_muts_pos]:
				max_muts_pos = pos
		#the pos of the max_muts_pos is declared as the reference for the distance of the remaining muts..
		score = 0
		for pos in range(cl_init, cl_end + 1):
			if pos in gene_accum_pos_dict:
				muts_perc_pos = float(gene_accum_pos_dict[pos])/ gene_mutations
				d = abs(pos - max_muts_pos)
				score += muts_perc_pos/ (CLUSTER_SCORE_CORRECTION **d)
		#once the cluster score is calculated, update the dict
		gene_clustering_score_dict[cluster_id] = score
	return gene_clustering_score_dict

'''
Compare scores of each gene (for non syn) with the distribution of scores observed among syn_mutations
'''
def get_gene_cluster_scores_z_dict(gene_cluster_scores_dict, gene_external_scores_dict):

	gene_cluster_scores_external_z_dict = {}
	#first i get all the values of the clusters of the external dataset
	scores_l = []

	for gene in gene_external_scores_dict:
		scores_l.append(gene_external_scores_dict[gene])

	scores_mean, scores_sd = calculate_mean_and_sd(scores_l)

	#then, i compare each value with the overall distribution
	for gene in gene_cluster_scores_dict:
		gene_cluster_scores_external_z_dict[gene] = get_z(gene_cluster_scores_dict[gene], scores_mean, scores_sd)


	return gene_cluster_scores_external_z_dict

'''
Given a dict{gene:{pos:num_accum_mutations}}
return the list of positions having more mutations than the number indicated by the threshold argument
'''
def get_significant_gene_positions(gene_accum_mut_pos_dict, minimum_mut_per_position_threshold):
	significant_gene_positions_l = []
	for pos in gene_accum_mut_pos_dict:
		if gene_accum_mut_pos_dict[pos] >=  minimum_mut_per_position_threshold:
			significant_gene_positions_l.append(pos) 
	return significant_gene_positions_l


def group_meaningful_positions(gene, gene_meaningful_positions):
	sorted_gene_meaningful_positions = sorted(gene_meaningful_positions)
	cluster_id = 1
	tmp_cluster_dict = {cluster_id: [sorted_gene_meaningful_positions[0]]}
	for i in range(1, len(sorted_gene_meaningful_positions)):
		if sorted_gene_meaningful_positions[i] <= sorted_gene_meaningful_positions[i-1] + INTRA_CLUSTER_MAX_DISTANCE:
			tmp_cluster_dict[cluster_id].append(sorted_gene_meaningful_positions[i])
		else:	
			cluster_id += 1
			tmp_cluster_dict[cluster_id] = [sorted_gene_meaningful_positions[i]]
	#dict = {cluster_id: (init_position, end_position)}
	gene_cluster_dict = {}
	for cluster_id in tmp_cluster_dict:
		gene_cluster_dict[cluster_id] = [tmp_cluster_dict[cluster_id][0], tmp_cluster_dict[cluster_id][len(tmp_cluster_dict[cluster_id])-1]]
	return gene_cluster_dict


def calculate_max_limit_cluster(gene, gene_meaningful_cluster_dict, cluster_id, pos, direction):
		if cluster_id == 1 and direction == 'left':
			return 0
		if cluster_id == len(gene_meaningful_cluster_dict.keys()) and direction == 'right':
			return int(cds_dict[gene])/3 if gene in cds_dict else pos
		elif direction == 'left':
			return pos - ( (pos - gene_meaningful_cluster_dict[cluster_id -1][1]) / 2)
		elif direction == 'right':
			return pos + ( (gene_meaningful_cluster_dict[cluster_id +1][0] - pos) / 2)

def calculate_limit_cluster(gene, pos, limit_pos, gene_mut_pos_l, direction):
	if direction == 'right':
		sorted_mut_pos = sorted(gene_mut_pos_l)
		for mut_pos in sorted_mut_pos:
			if mut_pos > limit_pos:
				return pos
			if pos < mut_pos <= pos + INTRA_CLUSTER_MAX_DISTANCE:
				pos = mut_pos
		return pos
	elif direction == 'left':
		reversed_mut_pos = sorted(gene_mut_pos_l, reverse = True)
		for mut_pos in reversed_mut_pos:
			if mut_pos < limit_pos:
				return pos
			if pos - INTRA_CLUSTER_MAX_DISTANCE <= mut_pos < pos:
				pos = mut_pos
		return pos

'''
Receive a dict with the meaningful positions already grouped ({cluster_id: [pos_init, pos_end]}) and it tries to extend
'''
def extend_meaningful_cluster(gene, gene_meaningful_cluster_dict, gene_accum_mut_pos_dict):
	gene_extended_meaningful_cluster_dict = {}

	for cluster_id in sorted(gene_meaningful_cluster_dict.keys()):
		init, end = gene_meaningful_cluster_dict[cluster_id][0], gene_meaningful_cluster_dict[cluster_id][1]

		limit_init = calculate_max_limit_cluster(gene, gene_meaningful_cluster_dict, cluster_id, init, 'left')
		limit_end = calculate_max_limit_cluster(gene, gene_meaningful_cluster_dict, cluster_id, end, 'right')

		init = calculate_limit_cluster(gene, init, limit_init, gene_accum_mut_pos_dict.keys(), 'left')
		end = calculate_limit_cluster(gene, end, limit_end, gene_accum_mut_pos_dict.keys(), 'right')

		gene_extended_meaningful_cluster_dict[cluster_id] = [init, end]

	return gene_extended_meaningful_cluster_dict


###################################################3
##
##   FUNCTION_CALLING
###################################################3

def run_oncoclust(options):

	syn_mut_fp, non_syn_mut_fp = options.syn_mut_fp, options.non_syn_mut_fp
	pos_pos, min_gene_mutations = options.pos_pos, options.min_gene_mutations
	output_fp = options.output_fp

	print '\nParsing the input files..'
	#dict = {gene: {pos: num_mutations}}
	non_syn_accum_mut_pos_dict = pos_f2accumulated_pos_dict(non_syn_mut_fp, pos_pos)
	syn_accum_mut_pos_dict = pos_f2accumulated_pos_dict(syn_mut_fp, pos_pos)

	# (A) Define the mutation clusters
	#dict = {gene:{cluster_id:[pos_init, pos_end]}}
	print '\nGrouping the lowly expected mutations into clusters..'
	print '\tNon_Synonimal mutations..'
	non_syn_cluster_coordinates_dict = cluster_mutations(non_syn_accum_mut_pos_dict, min_gene_mutations)
	print '\tCoding_silent mutations..'
	syn_cluster_coordinates_dict = cluster_mutations(syn_accum_mut_pos_dict, min_gene_mutations)

	#(B) Define the score of each cluster
	#-----> (B1) Calculate the mutations enclosed within each cluster
	#dict = {gene:{cluster_id: num_enclosed_muts}}
	print '\nGetting the number of mutations within each cluster...'
	print '\tNon_Synonimal mutations..'
	non_syn_cluster_muts_dict = get_enclosed_mutations(non_syn_cluster_coordinates_dict, non_syn_accum_mut_pos_dict)
	print '\tCoding_silent mutations..'
	syn_cluster_muts_dict = get_enclosed_mutations(syn_cluster_coordinates_dict, syn_accum_mut_pos_dict)

	#print 'NonSyn dict coordinates are: {0}'.format(non_syn_cluster_coordinates_dict['PIK3CA'])
	#print 'Enclosed mut dict is {0}'.format(non_syn_cluster_muts_dict['PIK3CA'])

	#-----> (B2) Calculate the score of each cluster: this score gives an idea of how specific are the mutation spatial concentration
	#dict = {gene :{cluster_id: clustering_score}}
	print '\nScoring the position specficity of the clusters..'
	print '\tNon_Synonimal mutations..'
	non_syn_cluster_scores_dict = score_clusters(non_syn_accum_mut_pos_dict, non_syn_cluster_coordinates_dict, non_syn_cluster_muts_dict)
	print '\tCoding_silent mutations..'
	syn_cluster_scores_dict = score_clusters(syn_accum_mut_pos_dict, syn_cluster_coordinates_dict, syn_cluster_muts_dict)

	#-----> Calculate the Z of each cluster score as compared to the scores of all the clusters of the CODING_SILENT muts
	#dict = {gene :{cluster_id: internal_z}}
	#print '\nRanking the NON_SYN cluster score per cluster by using results of the CODING_SILENT muts..'
	#non_syn_cluster_scores_external_z_dict = get_cluster_scores_external_z_dict(non_syn_cluster_scores_dict, syn_cluster_scores_dict)

	#-----> (B3) Retrieving the results per gene (i.e, combine the scores of all clusters of each gene)
	print '\nCombining info from cluster level to gene level..'
	non_syn_gene_cluster_scores_dict = group_clusters_of_genes(non_syn_cluster_scores_dict)
	syn_gene_cluster_scores_dict = group_clusters_of_genes(syn_cluster_scores_dict)

	# (C) Calculate the Z of the gene clustering score of non_syn vs syn (i.e. background model)
	print '\nObtaining the Z values by comparing gene cluster scores prot_affecting vs synonimous..'
	non_syn_gene_cluster_scores_z_dict = get_gene_cluster_scores_z_dict(non_syn_gene_cluster_scores_dict, syn_gene_cluster_scores_dict)

	# (D) Output
	print '\nOutputting the results..'
	create_output_file(options, non_syn_accum_mut_pos_dict, non_syn_cluster_coordinates_dict, \
	                   syn_gene_cluster_scores_dict, non_syn_gene_cluster_scores_dict, \
	                   non_syn_cluster_muts_dict, non_syn_gene_cluster_scores_z_dict)

	return



###################################################3
##
##   MAIN
###################################################3
'''
Note that I need as input two tdm files, in which each entry is a mutation.
I assume that the gene position is the first field, and I work with HUGO_symbols.
I assume that the AA_position of the mutation is the last field (this can be stated differently by user, see parser options)
'''
if __name__ == "__main__":

	#states the max number of aminoacids between 2 positions to group them into the same cluster
	INTRA_CLUSTER_MAX_DISTANCE = 5
	#the (minimum) probability threshold to consider the amount of mutations of a position as significant
	NOISE_BINOM_P_THRESHOLD = 0.01
	#the minimum mutations per position to consider that the position is meaningful
	MIN_MUTS_PER_MEANINGF_POS = 2
	#each position is weighted by correction^d, so the more widespread is the mut distribution within the cluster, the lower is the score
	CLUSTER_SCORE_CORRECTION = 1.4142


	#path to the cds fp (used for annotation purposes)
	working_folder = os.path.dirname(os.path.realpath(__file__))

	cds_fp = os.path.join(os.path.dirname(working_folder), 'data', 'max_len_transcript.tsv')
	cds_dict = f2dict(cds_fp, 1, -1)
	cancer_census_fp = os.path.join(os.path.dirname(working_folder), 'data','CGC_phenotype.tsv')

	options = parse_options()

	run_oncoclust(options)

	print '\n\nFinished!\n\n'


