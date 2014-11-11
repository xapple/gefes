#!/usr/bin/env python2

"""
A dirty script to compare two coverage matrices.
"""

# Modules #
import pandas

###############################################################################
theirs = pandas.io.parsers.read_csv("/Files/PHD/Server/Quince/homes/strainhack/data/StrainMock/coverage.tsv", sep='\t', index_col=0, encoding='utf-8')
ours = pandas.io.parsers.read_csv("/Files/PHD/Server/Quince/GEFES/views/projects/strain_mock/bins/concoct/coverage.tsv", sep='\t', index_col=0, encoding='utf-8')

theirs = theirs.loc[ours.index]
theirs.rename(columns=lambda x: x[16:], inplace=True)
theirs.values == ours.values