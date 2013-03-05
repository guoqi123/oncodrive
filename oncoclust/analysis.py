import logging

from utils import *

class OncoClustAnalysis(object):

	#states the max number of aminoacids between 2 positions to group them into the same cluster
	INTRA_CLUSTER_MAX_DISTANCE = 5

	#the (minimum) probability threshold to consider the amount of mutations of a position as significant
	#BINOMIAL_SIG_THRESHOLD_PER_POSITION = 0.01
	NOISE_BINOM_P_THRESHOLD = 0.01

	#each position is weighted by correction^d, so the more widespread is the mut distribution within the cluster, the lower is the score
	CLUSTER_SCORE_CORRECTION = 1.4142

	#the minimum mutations per position to consider that the position is meaningful
	MIN_MUTS_PER_MEANINGF_POS = 2


	def __init__(self):
		self.log = logging.getLogger("oncoclust.analysis")

	def get_gene_clusters(self, gene, gene_len, gene_muts, gene_accum_mut_pos_dict, cds_dict):
		'''
		First, it calculate the meaningful positions which are those with the minimum number of positions that
		are expected by less than 1% of probability according to the binomial model
		   -->	input : dict = {pos: num_mutations}  and gene_len, gene_muts
		   -->  output: dict = {cluster_id: [positions]}
		'''

		gene_cluster_dict = {}

		#The 'meaninful positions are those that have less than 1% of being found by chance. This is a binomial
		# cumulated probability model in which p(X >=n) = 1 - p(X <=(n-1)), where the parameters are: prob_of_success = 1/len_gene
		# and num_of_success = n, and num_of_trials = num_of_mutations
		#print '\nGetting gene meaningful positions..'

		#call func to return the minimum number of mutations per position to be included in the cluster
		# (note that 'self.NOISE_BINOM_P_THRESHOLD' is a global constant)
		minimum_mut_per_position_threshold = get_binomial_minimum_mut_per_position_threshold(gene_len, gene_muts, self.NOISE_BINOM_P_THRESHOLD, self.MIN_MUTS_PER_MEANINGF_POS)

		#given the minimum # of muts per position to be included in a cluster, retrieve the positions
		significant_gene_positions_l = self.get_significant_gene_positions(gene_accum_mut_pos_dict, minimum_mut_per_position_threshold)

		#if no meaningful positions, no cluster can be obtained
		if not len(significant_gene_positions_l) > 0:
			return {}

		#given the significant positions, group them
		# (applying the constraint of NOISE_BINOM_P_THRESHOLD, which is a global variable)
		# gene_dict = {cluster_id: [meaningful_pos_init, meaningful_pos_end]}
		gene_meaningful_cluster_dict = self.group_meaningful_positions(gene, significant_gene_positions_l)

		#given the meaningful positions, check if there is some non meaningful position around
		gene_extended_meaningful_cluster_dict = self.extend_meaningful_cluster(gene, gene_meaningful_cluster_dict, gene_accum_mut_pos_dict, cds_dict)

		#print 'Gene {0} and cluster dict is {1}'.format(gene, gene_meaningful_cluster_dict)
		#print 'Gene {0} and cluster extended dict is {1}'.format(gene, gene_extended_meaningful_cluster_dict)

		return gene_extended_meaningful_cluster_dict

	def pos_f2accumulated_pos_dict(self, mut_fp, pos_pos):
		'''
		returns a dict with the number (absolute) of mutations
		in each position of each gene: {gene: {pos:accumulated_mutations}}
		'''

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


	def cluster_mutations(self, accum_mut_pos_dict, min_gene_mutations, cds_dict):
		'''
		(A) Define meaningful positions in each gene, i.e. those positions showing a non-syn mutation occurrence higher than expected by chance.

		(B) The meaningful positions are clustered into groups that respect a maximal distance between two meaningful positions of the cluster

		Note that the input is a dict with the number of mutations per position: {gene: {pos:accumulated_mutations}}

		And returns a duct with the cluster_coordinates per gene (and per cluster): {gene: {cluster_id: [init_position, end_position]}}
		'''

		genes_sorted = sorted(accum_mut_pos_dict.keys())

		cluster_coordinates_dict = {}
		#---> Retrieval of clusters per gene
		for gene in genes_sorted:
			#---> (1) get the number of mutations of the gene
			gene_muts = sum([accum_mut_pos_dict[gene][pos] for pos in accum_mut_pos_dict[gene].keys()])

			# if there are few muts, as stated by the user (default =5),the gene is not included
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
			cluster_coordinates_dict[gene] = self.get_gene_clusters(gene, gene_len, gene_muts, accum_mut_pos_dict[gene], cds_dict)

		return cluster_coordinates_dict

	def get_enclosed_mutations(self, cluster_coordinates_dict, accum_mut_pos_dict):
		'''
		 inputs:
		  #dict = {gene: {cluster_id: [init_position, end_position]}}
		  #dict = {gene:{pos: num_mutations}}

		And returns: dict = {gene: {cluster_id: num_enclosed_mutations}}
		  (note that within the enclosed mutations, we take into account those that have been considered as being
		  in peak positions)
		'''

		cluster_summary_dict = {}
		for gene in cluster_coordinates_dict.keys():
			tmp_dict = {}
			for cluster_id in cluster_coordinates_dict[gene]:
				n_mutations = 0
				for position in accum_mut_pos_dict[gene]:
					cluster_init, cluster_end = cluster_coordinates_dict[gene][cluster_id][0], \
												cluster_coordinates_dict[gene][cluster_id][1]
					if position in range(cluster_init, cluster_end +1):
						n_mutations += accum_mut_pos_dict[gene][position]
				tmp_dict[cluster_id] = n_mutations
			cluster_summary_dict[gene] = tmp_dict
		return cluster_summary_dict


	def score_clusters(self, accum_mut_pos_dict, cluster_coordinates_dict, cluster_muts_dict):
		'''
		Calculate the scores for each cluster:
		- input: {gene:{cluster_id: [pos_init, pos_end]}}, {gene:{pos: muts}}
		- output: {gene: {cluster_id: score}}
		'''

		cluster_scores_dict = {}
		genes_sorted = sorted(cluster_coordinates_dict.keys())
		for gene in genes_sorted:
			#---> (1) get the number of mutations of the gene
			gene_muts = sum([accum_mut_pos_dict[gene][pos] for pos in accum_mut_pos_dict[gene].keys()])
			cluster_scores_dict[gene] = self.get_cluster_score_dict(gene, gene_muts,
															   accum_mut_pos_dict[gene],
															   cluster_coordinates_dict[gene])
		return cluster_scores_dict

	def group_clusters_of_genes(self, cluster_scores_dict):
		'''
		The clustering score of a gene is the sum of the scores of all the clusters of that gene
		'''

		gene_cluster_scores_dict = {}
		for gene in cluster_scores_dict:
			gene_cluster_scores_dict[gene] = sum(cluster_scores_dict[gene].values())
		return gene_cluster_scores_dict

	def get_cluster_score_dict(self, gene, gene_mutations, gene_accum_pos_dict, gene_cluster_coordinates_dict):
		'''
		Inputs:
			dict = {gene: {cluster_id: [pos_init, pos_end]}}
			dict = {gene: {pos: num_mutations}}
		The score here is calculated as the sum of the percentage of the mutations in each cluster pos, where each is divided
		by the distance from the pos with maximal mutations
		'''

		gene_clustering_score_dict = {}
		for cluster_id in gene_cluster_coordinates_dict:
			cl_init, cl_end = gene_cluster_coordinates_dict[cluster_id][0:2]

			#look for the pos of the cluster with more mutations
			max_muts_pos = cl_init
			for pos in range(cl_init, cl_end + 1):
				if pos in gene_accum_pos_dict and gene_accum_pos_dict[pos] > gene_accum_pos_dict[max_muts_pos]:
					max_muts_pos = pos
			#the pos of the max_muts_pos is declared as the reference for the distance of the remaining muts..
			score = 0
			for pos in range(cl_init, cl_end + 1):
				if pos in gene_accum_pos_dict:
					muts_perc_pos = float(gene_accum_pos_dict[pos]) / gene_mutations
					d = abs(pos - max_muts_pos)
					score += muts_perc_pos / (self.CLUSTER_SCORE_CORRECTION **d)
			#once the cluster score is calculated, update the dict
			gene_clustering_score_dict[cluster_id] = score
		return gene_clustering_score_dict

	def get_gene_cluster_scores_z_dict(self, gene_cluster_scores_dict, gene_external_scores_dict):
		'''
		Compare scores of each gene (for non syn) with the distribution of scores observed among syn_mutations
		'''

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

	def get_significant_gene_positions(self, gene_accum_mut_pos_dict, minimum_mut_per_position_threshold):
		'''
		Given a dict{gene:{pos:num_accum_mutations}}
		return the list of positions having more mutations than the number indicated by the threshold argument
		'''

		significant_gene_positions_l = []
		for pos in gene_accum_mut_pos_dict:
			if gene_accum_mut_pos_dict[pos] >=  minimum_mut_per_position_threshold:
				significant_gene_positions_l.append(pos)
		return significant_gene_positions_l

	def group_meaningful_positions(self, gene, gene_meaningful_positions):
		sorted_gene_meaningful_positions = sorted(gene_meaningful_positions)
		cluster_id = 1
		tmp_cluster_dict = {cluster_id: [sorted_gene_meaningful_positions[0]]}
		for i in range(1, len(sorted_gene_meaningful_positions)):
			if sorted_gene_meaningful_positions[i] <= sorted_gene_meaningful_positions[i-1] + self.INTRA_CLUSTER_MAX_DISTANCE:
				tmp_cluster_dict[cluster_id].append(sorted_gene_meaningful_positions[i])
			else:	
				cluster_id += 1
				tmp_cluster_dict[cluster_id] = [sorted_gene_meaningful_positions[i]]
		#dict = {cluster_id: (init_position, end_position)}
		gene_cluster_dict = {}
		for cluster_id in tmp_cluster_dict:
			gene_cluster_dict[cluster_id] = [tmp_cluster_dict[cluster_id][0], tmp_cluster_dict[cluster_id][len(tmp_cluster_dict[cluster_id])-1]]
		return gene_cluster_dict

	def calculate_max_limit_cluster(self, gene, gene_meaningful_cluster_dict, cluster_id, pos, direction, cds_dict):
			if cluster_id == 1 and direction == 'left':
				return 0
			if cluster_id == len(gene_meaningful_cluster_dict.keys()) and direction == 'right':
				return int(cds_dict[gene])/3 if gene in cds_dict else pos
			elif direction == 'left':
				return pos - ( (pos - gene_meaningful_cluster_dict[cluster_id -1][1]) / 2)
			elif direction == 'right':
				return pos + ( (gene_meaningful_cluster_dict[cluster_id +1][0] - pos) / 2)

	def calculate_limit_cluster(self, gene, pos, limit_pos, gene_mut_pos_l, direction):
		if direction == 'right':
			sorted_mut_pos = sorted(gene_mut_pos_l)
			for mut_pos in sorted_mut_pos:
				if mut_pos > limit_pos:
					return pos
				if pos < mut_pos <= pos + self.INTRA_CLUSTER_MAX_DISTANCE:
					pos = mut_pos
			return pos
		elif direction == 'left':
			reversed_mut_pos = sorted(gene_mut_pos_l, reverse = True)
			for mut_pos in reversed_mut_pos:
				if mut_pos < limit_pos:
					return pos
				if pos - self.INTRA_CLUSTER_MAX_DISTANCE <= mut_pos < pos:
					pos = mut_pos
			return pos

	'''
	Receive a dict with the meaningful positions already grouped ({cluster_id: [pos_init, pos_end]}) and it tries to extend
	'''
	def extend_meaningful_cluster(self, gene, gene_meaningful_cluster_dict, gene_accum_mut_pos_dict, cds_dict):
		gene_extended_meaningful_cluster_dict = {}

		for cluster_id in sorted(gene_meaningful_cluster_dict.keys()):
			init, end = gene_meaningful_cluster_dict[cluster_id][0], gene_meaningful_cluster_dict[cluster_id][1]

			limit_init = self.calculate_max_limit_cluster(gene, gene_meaningful_cluster_dict, cluster_id, init, 'left', cds_dict)
			limit_end = self.calculate_max_limit_cluster(gene, gene_meaningful_cluster_dict, cluster_id, end, 'right', cds_dict)

			init = self.calculate_limit_cluster(gene, init, limit_init, gene_accum_mut_pos_dict.keys(), 'left')
			end = self.calculate_limit_cluster(gene, end, limit_end, gene_accum_mut_pos_dict.keys(), 'right')

			gene_extended_meaningful_cluster_dict[cluster_id] = [init, end]

		return gene_extended_meaningful_cluster_dict

	###################################################3
	##
	##   FUNCTION_CALLING
	###################################################3

	def run(self, syn_mut_fp, non_syn_mut_fp, pos_pos, min_gene_mutations, cds_dict):

		self.log.info("Parsing the input files ...")

		#dict = {gene: {pos: num_mutations}}
		non_syn_accum_mut_pos_dict = self.pos_f2accumulated_pos_dict(non_syn_mut_fp, pos_pos)
		syn_accum_mut_pos_dict = self.pos_f2accumulated_pos_dict(syn_mut_fp, pos_pos)

		# (A) Define the mutation clusters
		#dict = {gene:{cluster_id:[pos_init, pos_end]}}
		self.log.info("Grouping the lowly expected mutations into clusters ...")

		self.log.info("  Non_Synonimal mutations ...")
		non_syn_cluster_coordinates_dict = self.cluster_mutations(non_syn_accum_mut_pos_dict, min_gene_mutations, cds_dict)

		self.log.info("  Coding_silent mutations ...")
		syn_cluster_coordinates_dict = self.cluster_mutations(syn_accum_mut_pos_dict, min_gene_mutations, cds_dict)

		#(B) Define the score of each cluster
		#-----> (B1) Calculate the mutations enclosed within each cluster
		#dict = {gene:{cluster_id: num_enclosed_muts}}
		self.log.info("Getting the number of mutations within each cluster ...")

		self.log.info("  Non_Synonimal mutations ...")
		non_syn_cluster_muts_dict = self.get_enclosed_mutations(non_syn_cluster_coordinates_dict, non_syn_accum_mut_pos_dict)

		self.log.info("  Coding_silent mutations ...")
		syn_cluster_muts_dict = self.get_enclosed_mutations(syn_cluster_coordinates_dict, syn_accum_mut_pos_dict)

		#-----> (B2) Calculate the score of each cluster: this score gives an idea of how specific are the mutation spatial concentration
		#dict = {gene :{cluster_id: clustering_score}}
		self.log.info("Scoring the position specficity of the clusters ...")

		self.log.info("  Non_Synonimal mutations ...")
		non_syn_cluster_scores_dict = self.score_clusters(non_syn_accum_mut_pos_dict, non_syn_cluster_coordinates_dict, non_syn_cluster_muts_dict)

		self.log.info("  Coding_silent mutations ...")
		syn_cluster_scores_dict = self.score_clusters(syn_accum_mut_pos_dict, syn_cluster_coordinates_dict, syn_cluster_muts_dict)

		#-----> Calculate the Z of each cluster score as compared to the scores of all the clusters of the CODING_SILENT muts
		#dict = {gene :{cluster_id: internal_z}}
		#print '\nRanking the NON_SYN cluster score per cluster by using results of the CODING_SILENT muts..'
		#non_syn_cluster_scores_external_z_dict = get_cluster_scores_external_z_dict(non_syn_cluster_scores_dict, syn_cluster_scores_dict)

		#-----> (B3) Retrieving the results per gene (i.e, combine the scores of all clusters of each gene)
		self.log.info("Combining info from cluster level to gene level ...")

		non_syn_gene_cluster_scores_dict = self.group_clusters_of_genes(non_syn_cluster_scores_dict)
		syn_gene_cluster_scores_dict = self.group_clusters_of_genes(syn_cluster_scores_dict)

		# (C) Calculate the Z of the gene clustering score of non_syn vs syn (i.e. background model)
		self.log.info("Obtaining the Z values by comparing gene cluster scores prot_affecting vs synonimous ...")
		non_syn_gene_cluster_scores_z_dict = self.get_gene_cluster_scores_z_dict(non_syn_gene_cluster_scores_dict, syn_gene_cluster_scores_dict)

		return (non_syn_accum_mut_pos_dict, non_syn_cluster_coordinates_dict,
			   syn_gene_cluster_scores_dict, non_syn_gene_cluster_scores_dict,
			   non_syn_cluster_muts_dict, non_syn_gene_cluster_scores_z_dict)
