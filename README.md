OncodriveCLUST
==============

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

Installation
------------

python setup.py install


