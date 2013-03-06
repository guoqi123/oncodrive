# OncodriveCLUST #

OncodriveCLUST is a method aimed to identify genes whose mutations are biased towards a large spatial clustering.
This method is designed to exploit the feature that mutations in cancer genes, especially oncogenes,
often cluster in particular positions of the protein. We consider this as a sign that mutations in these regions
change the function of these proteins in a manner that provides an adaptive advantage to cancer cells and
consequently are positively selected during clonal evolution of tumours, and this property can thus be used to nominate
novel candidate driver genes.

The method does not assume that the baseline mutation probability is homogeneous across all gene positions
but it creates a background model using silent mutations. Coding silent mutations are supposed to be under
no positive selection and may reflect the baseline clustering of somatic mutations. Given recent evidences
of non-random mutation processes along the genome, the assumption of homogenous mutation probabilities is likely
an oversimplication introducing bias in the detection of meaningful events.

## How it works ##

Detailed description is contained in the main manuscript. Briefly, the following steps are performed:
first, protein affecting mutations of each gene across a cohort of tumors are evaluated looking for
those protein residues having a number of mutations barely expected by chance. Second, these positions
are thereafter grouped to form mutation clusters. Third, each cluster is scored with a figure proportional
to the percentage of the gene mutations that are enclosed within that cluster and inversely related to its length.
The gene clustering score is obtained as the sum of the scores of all clusters (if any) found in that gene.
Finally, each gene clustering score is compared with the background model to obtain a significance value.
Background model is obtained performing the same steps than above but assessing only coding silent mutations.

## How it performs ##

We have analysed those entries of the COSMIC database annotated as whole gene screen as well as data provided
from 4 projects of the Cancer Genome Atlas. We demonstrated that the resulting candidate list of drivers is
strongly enriched by known cancer driver genes and particularly oncogenes, supporting the idea that this approach
can nominate novel driver candidates. In addition, comparison with methods based on other criteria
(namely, functional impact and mutation recurrence across the tumor cohort) demonstrated that the clustering
approach identifies known cancer drivers not detected by any of the other two methods, stressing the fact that
the combination of methods is beneficial to identify cancer drivers. We conclude that OncodriveCLUST is a method
that may be useful to identify cancer drivers through the assessment of the mutation clustering property that
may be complementary to other methods aimed to identify genes involved in the disease.

## Installation ##

OncodriveCLUST depends on some external libraries, [numpy](http://www.numpy.org/), [scipy](http://www.scipy.org/),
[pandas](http://pandas.pydata.org/) and [statsmodels](http://statsmodels.sourceforge.net/).

It is not easy to install them, specially scipy. We will give some clues on how to install them the way it works
for us but feel free to find your way. Once they are installed it is very easy to get OncodriveCLUST ready to work.

We recommend to use [virtualenv](http://www.virtualenv.org/). virtualenv is a tool to create isolated Python environments.
The basic problem being addressed is one of dependencies and versions, and indirectly permissions.
With virtualenv you can install the libraries and programs without having to be root.

If you are on *Mac OS X* or *Linux*, chances are that one of the following two commands will work for you:

	$ sudo easy_install virtualenv

or even better:

	$ sudo pip install virtualenv

One of these will probably install *virtualenv* on your system. Maybe it’s even in your package manager.
If you use *Ubuntu*, try:

	$ sudo apt-get install python-virtualenv

If you are on *Windows* and don’t have the *easy_install* command, you must install it first.
Check the *pip* and *distribute* on Windows section for more information about how to do that.
Once you have it installed, run the same commands as above, but without the sudo prefix.

Once you have virtualenv installed, just fire up a shell and create your own environment.

	$ virtualenv env

Now, whenever you want to work on a project, you only have to activate the corresponding environment.
On OS X and Linux, do the following:

	$ source env/bin/activate

If you are a Windows user, the following command is for you:

	$ env\scripts\activate

Either way, you should now be using your virtualenv (notice how the prompt of your shell has changed
to show the active environment).

Now you can just enter the following commands to get the OncodriveCLUST dependencies installed in your virtualenv:

	$ pip install -U distribute
	$ pip install -U numpy==1.6.1
	$ pip install -U scipy==0.9.0
	$ pip install -U pandas==0.10.1
	$ pip install -U statsmodels==0.4.0

Note that it would take quite long as they need to be compiled. One problem that could arise is that scipy require
BLAS and LAPACK or ATLAS libraries to be installed. In case they are not you have to download and compile them by yourself.
There is an installation guide at [http://www.scipy.org/Installing_SciPy](http://www.scipy.org/Installing_SciPy)

Then to get OncodriveCLUST installed run the following command:

	$ pip install https://bitbucket.org/bbglab/oncodriveclust/get/0.1.1.tar.gz

And that's all. The following command will allow you to check that is correctly installed by showing the command help:

	$ oncodriveclust --help
	-h, --help              show this help message and exit
	--version               show program's version number and exit
	-s PATH, --syn PATH     The path to the Synonimous mutations file to construct
	                        the background model
	-n PATH, --nonsyn PATH
	                        The path to the NON-Synonimous mutations file to be
							checked
	-p INT, --pos INT       AA position column index ('-1' by default, i.e. the last one)
	-m INT, --muts INT      Minimum number of mutations of a gene to be included
	                        in the analysis ('5' by default)
	-c, --coord             Use this argument for printing cluster coordinates in
	                        the output file
	-o PATH, --out PATH     Define the output file path

Note that -s, -n and -o are obligatory arguments

### Running an example ###

OncodriveCLUST requires two input files, one must contain the protein affecting mutations of the data set under analysis,
the other must contain the coding silent mutations that the method will use to construct the background model.
Both of them will be by default parsed as containing the gene symbol (i.e. HUGO) in the first column and the protein
position in the last column. Any other column between the first and the last will be ignored. Note that each entry
is assumed to be a different mutation (i.e, the mutation for a particular gene in a particular sample of the tumor
cohort).

With the [source code](https://bitbucket.org/bbglab/oncodriveclust/get/0.1.1.tar.gz) there is an example included.
We also provide with more [examples in a separate compressed file](https://bitbucket.org/bbglab/oncodriveclust/downloads/oncodriveclust-examples.tar.gz).
To run the default example:

	$ oncodriveclust -n examples/tcga.BRCA.nonsyn.txt -s examples/tcga.BRCA.syn.txt -o output_path -m 3

This will analyse the BRCA (breast invasive carcinoma) data set and the OncodriveCLUST output will be placed at
the output_path specified with the -o argument. Note that for the analyses of the manuscript, we have used an
arbitrary value of at least 10 mutations across the tumor cohort to include the gene in the analysis for the
COSMIC data set and a figure of 3 for the TCGA data sets analysis. This can be defined by using the -m argument.
With the optional flag -c the coordinates of the mutation clusters will be included in the output file.