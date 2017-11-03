#!/usr/bin/env python
import io
import os
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command

# Package meta-data.
NAME = 'RBNS_pipeline'
DESCRIPTION = 'A pipeline to process RBNS data'
URL = 'https://bitbucket.org/pfreese/rbns_pipeline'
EMAIL = 'freese.peter@gmail.com'
AUTHOR = 'Peter Freese'

# What packages are required for this module to be executed?
REQUIRED = [
            # 'requests', 'maya', 'records',
            'setuptools',
            'matplotlib',
            'Pillow',
            'PyPDF2'
            'forgi',
            'numpy',
            'pdfnup',
            'scipy',
            'simplejson']

here = os.path.abspath(os.path.dirname(__file__))

# Where the magic happens:
setup(
    name=NAME,
    version='1.0',
    description=DESCRIPTION,
    # long_description=long_description,
    author=AUTHOR,
    author_email=EMAIL,
    url=URL,
    packages=find_packages(exclude=('tests',)),
                install_requires=REQUIRED,
                include_package_data=True,
                license='MIT',
                classifiers=[
                'License :: OSI Approved :: MIT License',
                'Programming Language :: Python',
                'Programming Language :: Python :: 2.7'] )






