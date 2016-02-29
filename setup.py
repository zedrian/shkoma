#!/usr/bin/env python
from setuptools import setup


description = "More information at http://github.com/zedrian/shkoma"


setup(
    name='shkoma',
    version='0.1',
    author='Mihail Shapiro, Artsiom Kopats',
    author_email='artsiom.kopats@gmail.com',
    url='http://github.com/zedrian/shkoma',
    description='Protein data analysis',
    long_description=description,
    license='MIT',
    install_requires=['uniprot', 'numpy', 'biopython'],
    packages=['shkoma']
)
