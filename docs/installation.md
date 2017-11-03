![Logo](../img/SHAPEware-logo.png)

***
# How to install RBNS pipeline dependencies

## Python

It is highly recommended that you use python through [Conda](https://conda.io/docs/user-guide/install/index.html).

- The simplejson package is also required ('pip install simplejson').

## Forgi

- Note that the forgi package has its own dependencies, as listed on [the forgi website](https://viennarna.github.io/forgi/download.html). On test installion, this required installing future ('pip install future') ad logging_exceptions ('pip install logging_exceptions').

## Weblogo

The easiest way to obtain [Weblogo](http://weblogo.threeplusone.com/manual.html) is via pip:

pip install weblogo

Weblogo must be able to be called from the command line. 'pip install weblogo' will automatically add it to your conda bin, but if you install it stand-alone, be sure to add weblogo executable to your $PATH. 


## RNAfold

RNAfold program can be downloaded from [Vienna RNA](https://www.tbi.univie.ac.at/RNA/#download). Once downloaded and installed, be sure to add the it to your $PATH (i.e. /path-to-RNAfold/bin should be on your path, so 'RNAfold' from anywhere on the command line begins the program).

