# OncodriveCLUST

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

## How it works

Detailed description is contained in the main manuscript. Briefly, the following steps are performed:
first, protein affecting mutations of each gene across a cohort of tumors are evaluated looking for
those protein residues having a number of mutations barely expected by chance. Second, these positions
are thereafter grouped to form mutation clusters. Third, each cluster is scored with a figure proportional
to the percentage of the gene mutations that are enclosed within that cluster and inversely related to its length.
The gene clustering score is obtained as the sum of the scores of all clusters (if any) found in that gene.
Finally, each gene clustering score is compared with the background model to obtain a significance value.
Background model is obtained performing the same steps than above but assessing only coding silent mutations.

## How it performs

We have analysed those entries of the COSMIC database annotated as whole gene screen as well as data provided
from 4 projects of the Cancer Genome Atlas. We demonstrated that the resulting candidate list of drivers is
strongly enriched by known cancer driver genes and particularly oncogenes, supporting the idea that this approach
can nominate novel driver candidates. In addition, comparison with methods based on other criteria
(namely, functional impact and mutation recurrence across the tumor cohort) demonstrated that the clustering
approach identifies known cancer drivers not detected by any of the other two methods, stressing the fact that
the combination of methods is beneficial to identify cancer drivers. We conclude that OncodriveCLUST is a method
that may be useful to identify cancer drivers through the assessment of the mutation clustering property that
may be complementary to other methods aimed to identify genes involved in the disease.

## Installation

OncodriveCLUST depends on some external libraries, [numpy](http://www.numpy.org/), [scipy](http://www.scipy.org/)
and [statsmodels](http://statsmodels.sourceforge.net/).

There are many ways to get everything ready to use OncodriveCLUST, but the most recommended method is using
[virtualenv](http://www.virtualenv.org/), so we will start explaining it. We will also cover how to install it
the old fashion way (*system-wide*) and will cover some possible problems that may arise with dependencies installation.

### virtualenv

Virtualenv is probably what you want to use during development, and if you have shell access to your production machines,
you will probably want to use it there, too.

What problem does virtualenv solve? When working with many Python programs the chance to have conflicts
between libraries required by different programs or the required version of Python increases.
Virtualenv enables multiple side-by-side installations of Python, one for each project. It doesn’t actually install
separate copies of Python, but it does provide a clever way to keep different project environments isolated.

If you are on *Mac OS X* or *Linux*, chances are that one of the following two commands will work for you:

> $ sudo easy_install virtualenv

or even better:

> $ sudo pip install virtualenv

One of these will probably install *virtualenv* on your system. Maybe it’s even in your package manager.
If you use *Ubuntu*, try:

> $ sudo apt-get install python-virtualenv

If you are on *Windows* and don’t have the *easy_install* command, you must install it first.
Check the *pip* and distribute on Windows section for more information about how to do that.
Once you have it installed, run the same commands as above, but without the sudo prefix.

Once you have virtualenv installed, just fire up a shell and create your own environment.

> $ virtualenv venv

Now, whenever you want to work on a project, you only have to activate the corresponding environment.
On OS X and Linux, do the following:

> $ source venv/bin/activate

If you are a Windows user, the following command is for you:

> $ venv\scripts\activate

Either way, you should now be using your virtualenv (notice how the prompt of your shell has changed
to show the active environment).

Now you can just enter the following command to get OncodriveCLUST installed in your virtualenv:

> $ pip install https://bitbucket.org/bbglab/oncodriveclust/get/0.2.tar.gz

That's all.

### System-wide

May be you prefer to do the thigs more manually or you found some problem following the previous steps.

Just download https://bitbucket.org/bbglab/oncodriveclust/get/0.2.tar.gz and uncompress it.
Go inside the uncompressed folder and run:

> $ python setup.py install

This should install all the required dependencies and install the tool. If there is any problem installing dependencies
then you should install them manually by your own and then try again. It is recommended to use the system package management system
(i.e. apt-get in Ubuntu or yum in fedora).

### Running the example

With the [source code](https://bitbucket.org/bbglab/oncodriveclust/get/0.2.tar.gz) there is included an example.
Download it and execute:

> $ oncodriveclust -n examples/tcga.BRCA.nonsyn.txt -s examples/tcga.BRCA.syn.txt -o output_path -m 3