#!/usr/bin/env python
from setuptools import setup


description = "More information at http://github.com/zedrian/shkoma"


setup(
    name='shkoma',
    version='0.1',
    author='Artsiom Kopats, Mihail Shapiro',
    author_email='artsiom.kopats@gmail.com, mshapira2016@gmail.com',
    url='http://github.com/zedrian/shkoma',
    description='Protein data analysis',
    long_description=description,
    license='MIT',
    install_requires=['uniprot', 'numpy', 'biopython', 'rpy2'],
    packages=['shkoma']
)
