#!/usr/bin/env python2

"""
A dirty script to compare two coverage matrices.
"""

# Modules #
import pandas

###############################################################################
df1 = pandas.io.parsers.read_csv("/Files/PHD/Server/Quince/homes/strainhack/data/StrainMock/coverage.tsv", sep='\t', index_col=0, encoding='utf-8')
df2 = pandas.io.parsers.read_csv("/Files/PHD/Server/Quince/GEFES/views/projects/strain_mock/bins/concoct/coverage.tsv", sep='\t', index_col=0, encoding='utf-8')
df0 = df1.loc[df2.index]
df0.rename(columns=lambda x: x[16:], inplace=True)
df0.values == df2.values