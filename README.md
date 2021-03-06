# OncodriveCLUST #

## --> Warning: OncodriveCLUST is outdated! ##

We have developed a new clustering algorithm, OncodriveCLUSTL, that outperforms OncodriveCLUST in the detection of protein-coding cancer driver genes and can be applied to non-coding regions of the genome. You can read OncodriveCLUSTL publication with the complete benchmark and details: Arnedo-Pac C, et al. OncodriveCLUSTL: a sequence-based clustering method to identify cancer drivers. Bioinformatics. 2019;35(22):4788–4790. [doi:10.1093/bioinformatics/btz501](https://doi.org/10.1093/bioinformatics/btz501)

OncodriveCLUSTL is available as an installable Python 3.5 package at pip and conda. The source code and running examples are freely available at [bitbucket.org/bbglab/oncodriveclustl](https://bitbucket.org/bbglab/oncodriveclustl) under GNU Affero General Public License. There is a web version at [bbglab.irbbarcelona.org/oncodriveclustl](http://bbglab.irbbarcelona.org/oncodriveclustl/home)

## Description ##

OncodriveCLUST is a method aimed to identify genes whose mutations are biased towards a large spatial clustering. This method is designed to exploit the feature that mutations in cancer genes, especially oncogenes, often cluster in particular positions of the protein. We consider this as a sign that mutations in these regions change the function of these proteins in a manner that provides an adaptive advantage to cancer cells and consequently are positively selected during clonal evolution of tumours, and this property can thus be used to nominate novel candidate driver genes.

The method does not assume that the baseline mutation probability is homogeneous across all gene positions but it creates a background model using silent mutations. Coding silent mutations are supposed to be under no positive selection and may reflect the baseline clustering of somatic mutations. Given recent evidences of non-random mutation processes along the genome, the assumption of homogenous mutation probabilities is likely an oversimplication introducing bias in the detection of meaningful events.

## How it works ##

Detailed description is contained in the main manuscript. Briefly, the following steps are performed: first, protein affecting mutations of each gene across a cohort of tumors are evaluated looking for those protein residues having a number of mutations barely expected by chance. Second, these positions are thereafter grouped to form mutation clusters. Third, each cluster is scored with a figure proportional to the percentage of the gene mutations that are enclosed within that cluster and inversely related to its length. The gene clustering score is obtained as the sum of the scores of all clusters (if any) found in that gene. Finally, each gene clustering score is compared with the background model to obtain a significance value. Background model is obtained performing the same steps than above but assessing only coding silent mutations.

## How it performs ##

We have analysed those entries of the COSMIC database annotated as whole gene screen as well as data provided from 4 projects of the Cancer Genome Atlas. We demonstrated that the resulting candidate list of drivers is strongly enriched by known cancer driver genes and particularly oncogenes, supporting the idea that this approach can nominate novel driver candidates. In addition, comparison with methods based on other criteria (namely, functional impact and mutation recurrence across the tumor cohort) demonstrated that the clustering approach identifies known cancer drivers not detected by any of the other two methods, stressing the fact that the combination of methods is beneficial to identify cancer drivers. We conclude that OncodriveCLUST is a method that may be useful to identify cancer drivers through the assessment of the mutation clustering property that may be complementary to other methods aimed to identify genes involved in the disease.

## Installation ##

OncodriveCLUST depends on Python 3 and some external libraries, [numpy](http://www.numpy.org/), [scipy](http://www.scipy.org/), [pandas](http://pandas.pydata.org/) and [statsmodels](http://statsmodels.sourceforge.net/).

The easiest way to install all this software stack is using the well known [Anaconda Python distribution](http://continuum.io/downloads#py34).

Then to get OncodriveCLUST installed run the following command:

	(env) $ pip install oncodriveclust

And that's all. The following command will allow you to check that is correctly installed by showing the command help:

	(env) $ oncodriveclust --help

	usage: oncodriveclust [-h] [--version] [-o PATH] [--cgc PATH] [-m INT] [-c]
						  [-p INT]
						  NON-SYN-PATH SYN-PATH GENE-TRANSCRIPTS

	Run OncodriveCLUST analysis

	positional arguments:
	  NON-SYN-PATH         The path to the NON-Synonymous mutations file to be
						   checked
	  SYN-PATH             The path to the Synonymous mutations file to construct
						   the background model
	  GENE-TRANSCRIPTS     The path of a file containing transcripts length for
						   genes

    optional arguments:
      -h, --help            show this help message and exit
      --version             show program's version number and exit
      -o PATH, --out PATH   Define the output file path
      --cgc PATH            The path of a file containing CGC data
      -m INT, --muts INT    Minimum number of mutations of a gene to be included
                            in the analysis ('5' by default)
      -c, --coord           Use this argument for printing cluster coordinates in
                            the output file
      --pos INT             AA position column index ('-1' by default)
      -d INT, --dist INT    Intra cluster maximum distance ('5' by default)
      -p FLOAT, --prob FLOAT
                            Probability of the binomial model to find cluster
                            seeds ('0.01' by default)
      --dom PATH            The path of a file containing gene domains
      -L LEVEL, --log-level LEVEL
                            Define the loggging level


### Running an example ###

OncodriveCLUST requires three input files, one must contain the protein affecting mutations of the data set under analysis, the other must contain the coding silent mutations that the method will use to construct the background model and the last one contain the CDS length for the genes transcripts.

The first two files will be by default parsed as containing the gene symbol (i.e. *HUGO*) in the first column and the protein position in the last column. Any other column between the first and the last will be ignored. Note that each entry is assumed to be a different mutation (i.e, the mutation for a particular gene in a particular sample of the tumor cohort). The file containing CDS lengths have three columns: the gene symbol, the transcript symbol and the CDS length. OncodriveCLUST will get the gene CDS length from the longest transcript automatically.

With the [source code](https://bitbucket.org/bbglab/oncodriveclust/get/0.3.tar.gz) you will get some example files, a file containing transcripts CDS lengths and a file containing [CGC](http://www.sanger.ac.uk/genetics/CGP/Census/) phenotypes for known genes, so it is recommended to download and uncompress it (just click on the previous link, uncompress it and then change the current directory).

```
(env) $ curl -o oncodriveclust.tar.gz https://bitbucket.org/bbglab/oncodriveclust/get/0.3.tar.gz
(env) $ tar -xvf oncodriveclust.tar.gz
(env) $ cd bbglab-oncodriveclust-<code>
```

We also provide with more [examples in a separate compressed file](https://bitbucket.org/bbglab/oncodriveclust/downloads/oncodriveclust_examples.tar.gz).

To run the default example:

	(env) $ oncodriveclust -m 3 --cgc data/CGC_phenotype.tsv examples/tcga.BRCA.nonsyn.txt examples/tcga.BRCA.syn.txt data/gene_transcripts.tsv

This will analyse the BRCA (breast invasive carcinoma) data set and the OncodriveCLUST output will be placed at the path specified with the -o argument (if not specified then the current directory will be used). Note that for the analyses of the manuscript, we have used an arbitrary value of at least 10 mutations across the tumor cohort to include the gene in the analysis for the COSMIC data set and a figure of 3 for the TCGA data sets analysis. This can be defined by using the -m argument. With the optional flag -c the coordinates of the mutation clusters will be included in the output file.

The use of "--cgc data/CGC_phenotype.tsv" is optional, but if specified then the results will contain information about the CGC phenotype for each gene if available.

The results file will look something like:

	(env) $ head oncodriveclust-results.tsv
	GENE	CGC	GENE_LEN	GENE_NUM_MUTS	MUTS_IN_CLUST	NUM_CLUSTERS	GENE_SCORE	ZSCORE	PVALUE	QVALUE
    KRAS	Dom	190	3	3	1	1.0	4.17671455905	1.47874885544e-05	0.0005077037737
    KLRG2		410	3	3	1	1.0	4.17671455905	1.47874885544e-05	0.0005077037737
    KRT38		457	3	3	1	1.0	4.17671455905	1.47874885544e-05	0.0005077037737
    AKT1	Dom	481	12	11	1	0.916666666667	3.79776654893	7.30028674072e-05	0.00187982383574
    PIK3CA	Dom	1069	186	178	7	0.855624278036	3.52018384847	0.00021562388305	0.00444185199083
    VASN		674	5	4	1	0.8	3.26723933475	0.000543009128935	0.00932165671338
    KIF19		999	3	3	1	0.784521187523	3.1968513525	0.000694682523388	0.0100330867162
    ZNF598		905	4	3	1	0.75	3.03987052868	0.00118339934876	0.0100330867162
    FGFR2	Dom	823	4	3	1	0.75	3.03987052868	0.00118339934876	0.0100330867162

In addition, the user can add a column with the details of the clusters identified by OncodriveCLUST (i.e. aa_start and aa_end delimiting such clusters and the number of mutations enclosed in each) by using the -c argument. Finally, the domains of the protein encoded by each gene can be added as another column in the output by using the --dom argument: "--dom data/pfam_domains.txt". In this file, we provide the data retrieved from Pfam; the user can use another self-prepared file if desired by indicating the corresponding path to this file.