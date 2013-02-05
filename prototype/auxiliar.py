import math
import rpy2.robjects as robjects
from operator import itemgetter

###################################################3
##
##   Auxiliar 
###################################################3
'''
convert a file with a key-value to a dict
'''
def f2dict(f_p, k_position, v_position):
	dict = {}
	f = open(f_p, 'r')
	h_s = f.readline().rstrip().split("\t")
	for l in f:
		try:
			l_s = l.rstrip().split("\t")
			key, value = l_s[k_position], l_s[v_position]
			dict[key] = value
		except:
			continue
	f.close()
	return dict

'''
Given a len of gene and a number of gene mutations, I find the minimum number of mutations having a prob of occurring
in a position <= 1%

I assume that the probability of having n muts in a single position is like having n succceses in gene_muts trials, 
  where the prob_success is 1/gene_len
'''
def get_binomial_minimum_mut_per_position_threshold(gene_len, gene_muts, sig_cutoff):
	prob_success = 1/float(gene_len)
	num_trials = int(gene_muts)

	pbinom = robjects.r['pbinom']
	minimum_mut_per_position_threshold = 2

	for num_success in range(2, gene_muts):
		binom_p = pbinom(num_success, num_trials, prob_success)[0]

		if 1-binom_p <= sig_cutoff:
			minimum_mut_per_position_threshold = num_success
			break
	return minimum_mut_per_position_threshold


def get_z(value, mean, sd):
	return (float(value) - mean) / sd

'''
add a column with P value and another with FDR values to
the output_matrix (already ordered by Z)
'''
def add_fdr(matrix, pos= -1, sense = 'OVER'):

	pnorm = robjects.r['pnorm']
	padj = robjects.r['p.adjust']

	z = []
	[z.append(float(l[-1])) for l in matrix] if sense == 'UNDER' else [z.append(float(l[-1]) * -1) for l in matrix]

	rz = robjects.FloatVector(z)
	
	rp = pnorm(rz)
	rf = padj(rp, method = 'fdr')

	#add the calculated p and fdr values to the output matrix	
	for i in range(0, len(matrix)):
		matrix[i] += [rp[i], rf[i]]

	return matrix


def order_matrix(matrix, position, order = 'DEC'):
	if order == 'DEC':
		matrix.sort(key=itemgetter(position), reverse = True)
	else:
		matrix.sort(key=itemgetter(position), reverse = False)
	return matrix

def calculate_mean_and_sd(values_list):
	final_values_list = []
	for value in values_list:	
		try:
			final_values_list.append(float(value))
		except:
			continue
	n = len(final_values_list)
	mean = sum(final_values_list)/n
	sd = math.sqrt(sum((x-mean)**2 for x in final_values_list) / n)
	return mean, sd
