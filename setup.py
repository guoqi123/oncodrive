"""
OncoCLUST
=========

"""

import distribute_setup
distribute_setup.use_setuptools()

from setuptools import setup, find_packages

from oncoclust import VERSION, AUTHORS, AUTHORS_EMAIL

setup(
	name = "oncoclust",
	version = VERSION,
	packages = ["oncoclust"],
	scripts = [],

	install_requires = [
		"numpy>=1.6.1",
		"scipy>=0.11.0",
		"pandas>=0.10.1",
		"statsmodels>=0.4",
	],

	package_data = {
		# If any package contains *.txt or *.pdf files, include them:
		"" : ["*.txt", "*.md", "*.pdf", "*.html"]
	},

	entry_points = {
		'console_scripts': [
			'oncoclust = oncoclust.command:main'
		]
	},

	# metadata for upload to PyPI
	author = AUTHORS,
	author_email = AUTHORS_EMAIL,
	description = "OncoCLUST",
	license = "AGPL 3.0",
	keywords = "",
	url = "https://bitbucket.org/bbglab/oncoclust",
	long_description = __doc__,

	classifiers = [
		"Development Status :: 4 - Beta",
		"Intended Audience :: Bioinformatics",
		"License :: OSI Approved :: GNU Affero General Public License (AGPL)",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"Natural Language :: English",
		"Operating System :: OS Independent",
		"Programming Language :: Python",
		"Programming Language :: Python :: 2.7",
		"Topic :: Scientific/Engineering",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
