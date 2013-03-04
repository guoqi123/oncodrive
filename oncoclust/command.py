import os
import argparse
import logging

from . import VERSION
from analysis import OncoClustAnalysis
from utils import *

_log_level_map = {
	"debug" : logging.DEBUG,
	"info" : logging.INFO,
	"warn" : logging.WARN,
	"error" : logging.ERROR,
	"critical" : logging.CRITICAL,
	"notset" : logging.NOTSET }

class Command(object):
	def __init__(self, prog=None, desc=""):
		parser = argparse.ArgumentParser(prog=prog, description=desc)

		self._add_arguments(parser)

		parser.add_argument("-L", "--log-level", dest="log_level", metavar="LEVEL", default=None,
							choices=["debug", "info", "warn", "error", "critical", "notset"],
							help="Define the loggging level")

		self.args = parser.parse_args()

		logging.basicConfig(
			format = "%(asctime)s %(name)s %(levelname)-5s : %(message)s",
			datefmt = "%Y-%m-%d %H:%M:%S")

		if self.args.log_level is None:
			self.args.log_level = "info"
		else:
			self.args.log_level = self.args.log_level.lower()

		logging.getLogger("oncoclust").setLevel(_log_level_map[self.args.log_level])

		self.log = logging.getLogger("oncoclust")

	def _add_arguments(self, parser):
		pass

	def _check_args(self):
		pass

	def run(self):

		# Check arguments

		self._check_args()

class OncoClustCommand(Command):
	def __init__(self):
		Command.__init__(self, prog="oncoclust", desc="Run OncoCLUST")

		self.root_path = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))

		self.cds_fp = os.path.join(self.root_path, 'data', 'max_len_transcript.tsv')

		self.cancer_census_fp = os.path.join(self.root_path, 'data','CGC_phenotype.tsv')

	def _add_arguments(self, parser):
		Command._add_arguments(self, parser)

		parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

		parser.add_argument("-s", "--syn", dest="syn_mut_fp", metavar="PATH",
				help="The path to the Synonimous mutations file to construct the background model")

		parser.add_argument("-n", "--nonsyn", dest="non_syn_mut_fp", metavar="PATH",
				help="The path to the NON-Synonimous mutations file to be checked")
	
		parser.add_argument("-p", "--pos", dest="pos_pos", type=int, default=-1, metavar="INT",
				help="AA position column index ('-1' by default)")
	
		parser.add_argument("-m", "--muts", dest="min_gene_mutations", type=int, default=5, metavar="INT",
				help="Minimum number of mutations of a gene to be included in the analysis ('5' by default)")
	
		parser.add_argument("-o", "--out", dest="output_fp", metavar="PATH",
				help="Define the output file path")
	
		self.args = parser.parse_args()

	def _check_args(self):
		if not self.args.syn_mut_fp or not self.args.non_syn_mut_fp or\
				not os.path.exists(self.args.syn_mut_fp) or not os.path.exists(self.args.non_syn_mut_fp):
			self.log.error("Please specify a valid path for both Syn and Nonsyn mutations files by using the corresponding"
					 " arguments. Use --help for further details.")
			exit(-1)

		if not self.args.output_fp:
			self.log.error("Please specify a path for output file by using the corresponding"
					 " arguments. Use --help for further details.")
			exit(-1)

	def create_output_file(self, cds_dict, cgc_dict,
						   non_syn_accum_mut_pos_dict, non_syn_cluster_coordinates_dict,
						   syn_gene_cluster_scores_dict, non_syn_gene_cluster_scores_dict,
						   non_syn_cluster_muts_dict, non_syn_gene_cluster_scores_external_z_dict):

		#Output file
		f_out_fp = "".join([self.args.output_fp,
							 '.scoreCorrection' + str(OncoClustAnalysis.CLUSTER_SCORE_CORRECTION),
							 '.minMuts' + str(self.args.min_gene_mutations),
							 '.GENE.LEVEL.muts_clustering.txt'])

		self.log.debug("> {0}".format(f_out_fp))

		f_out = open(f_out_fp, 'w')
		f_out.write('\t'.join(['Gene', 'CGC', 'Gene_len', 'Gene_Muts', 'n_clusters', 'Muts_in_clusters', 'Gene_cluster_score','Z_value', 'P_value', 'Q_value' ]))

		m_out = []
		m_not_included_out = []
		for gene in non_syn_gene_cluster_scores_dict:
			cgc = cgc_dict[gene] if gene in cgc_dict else ''
			gene_len = int(cds_dict[gene])
			gene_muts = sum([non_syn_accum_mut_pos_dict[gene][pos] for pos in non_syn_accum_mut_pos_dict[gene].keys()])
			n_clusters = len(non_syn_cluster_coordinates_dict[gene].keys())

			if n_clusters> 0:
				muts_clusters = sum(non_syn_cluster_muts_dict[gene].values())
				gene_additive_cluster_score = non_syn_gene_cluster_scores_dict[gene]
				z_ext = non_syn_gene_cluster_scores_external_z_dict[gene] if gene in non_syn_gene_cluster_scores_external_z_dict else 'NA'
				m_out.append([gene, cgc, gene_len, gene_muts, n_clusters, muts_clusters, gene_additive_cluster_score, z_ext])
			else:
				m_not_included_out.append([gene, cgc, gene_len, gene_muts, 0, 0, 0, 'NA', 'NA', 'NA'])

		m_out = order_matrix(m_out, -1)
		m_out = add_fdr(m_out, -1)
		m_out += m_not_included_out

		for l_out in m_out:
			f_out.write('\n' + '\t'.join(str(e) for e in l_out))
		f_out.close()

	def run(self):
		'''
		Note that I need as input two tdm files, in which each entry is a mutation.
		I assume that the gene position is the first field, and I work with HUGO_symbols.
		I assume that the AA_position of the mutation is the last field (this can be stated differently by user, see parser options)
		'''

		Command.run(self)

		analysis = OncoClustAnalysis()

		self.log.info("Loading CDS data ...")

		#dict = {gene: length}
		cds_dict = f2dict(self.cds_fp, 1, -1)

		self.log.info("Running analysis ...")

		results = analysis.run(
			self.args.syn_mut_fp, self.args.non_syn_mut_fp,
			self.args.pos_pos, self.args.min_gene_mutations,
			cds_dict)

		self.log.info("Loading CGC data ...")

		cgc_dict = f2dict(self.cancer_census_fp, 0, 1)

		self.log.info("Saving results ...")

		self.create_output_file(cds_dict, cgc_dict, *results)


def main():
	"""
	Main entry point
	"""
	OncoClustCommand().run()

if __name__ == "__main__":
	main()