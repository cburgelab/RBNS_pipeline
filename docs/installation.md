![Logo](../img/RBNS_logo.png)

***
# How to install RBNS pipeline dependencies

## Python

It is highly recommended that you use python and other dependencies through a [Conda](https://conda.io/docs/user-guide/install/index.html) installation:

	git clone https://pfreese@bitbucket.org/pfreese/rbns_pipeline.git
	cd rbns_pipeline
	./install.sh

The install.sh should take a while while Conda installs the dependencies. After activating your Conda rbns_pipeline environment via: 'source activate rbns_pipeline', you can then install the remaining dependencies:

- The simplejson package is required ('pip install simplejson').

- The forgi package is required ('pip install forgi').
	- Note that the forgi package has its own dependencies, as listed on [the forgi website](https://viennarna.github.io/forgi/download.html). On test installion, this required installing future ('pip install future') and logging_exceptions ('pip install logging_exceptions').


- The pdfnup command line program is required ('pip install pdfnup'). 
	- Note that pdfnup.py script (which should be within your site-packages at /path/to/your/miniconda2/envs/rbns_pipeline/lib/python2.7/site-packages or /path/to/your/anaconda2/lib/python2.7/site-packages), may need to be modified in the following way: pdfnup.py has the line 'from pyPdf.pdf import PageObject, ImmutableSet, ContentStream', but the ImmutableSet import my throw an error; if that line is changed to just 'from pyPdf.pdf import PageObject, ContentStream' and then on a separate line add 'from sets import ImmutableSet', and everything should be OK. Just be sure to test via 'python pdfnup.py' and make sure there are no errors.

## Weblogo

The easiest way to obtain [Weblogo](http://weblogo.threeplusone.com/manual.html) is via pip:

	pip install weblogo

Weblogo must be able to be called from the command line. 'pip install weblogo' will automatically add it to your conda bin, but if you install it stand-alone, be sure to add the weblogo executable to your $PATH. 


## RNAfold

The RNAfold program should be installed via the above as the 'viennarna=2.3.5' Conda dependency. Check to make sure 'RNAfold' is on your path and can be called from anywhere the command line; if it isnot, it can be downloaded from [Vienna RNA](https://www.tbi.univie.ac.at/RNA/#download). Once downloaded and installed, be sure to add it to your $PATH (i.e. /path-to-RNAfold/bin should be on your path).

