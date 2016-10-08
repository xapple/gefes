#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script that reads an excel file and produces the correctly formated JSON files.

You can use this script from the shell like this:
$ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/micans_v6_exp1/metadata.xlsx
$ ~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/micans_v6_exp1_plexed/metadata_plexed.xlsx

NB: The excel sheet must have an non-empty cell in the position (0,0)
NB: The excel sheet must have an non-empty cell in the position (0,1)
"""

# Modules #
import sys, os, pandas, codecs, pystache

# Internal modules #
from correspondence import corr

# Third party modules #
import numpy

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class JsonTemplate(object):
    def __init__(self, content):
        self.content = content

###############################################################################
# Get the shell arguments #
if len(sys.argv) < 2: sys.exit(sys.modules[__name__].__doc__)
xlsx_path = sys.argv[1]

# Check that the path is valid #
if not os.path.exists(xlsx_path): raise Exception("No file at %s." % xlsx_path)

# Load data #
xlsx = pandas.ExcelFile(xlsx_path)
df = xlsx.parse(xlsx.sheet_names[0])

# Use the first row (of the dataframe) as the column index and then delete the first row #
df.columns = df.iloc[0]
df = df.drop(0)

# Main loop - One JSON file per sample #
for i, row in df.iterrows():

    # Using the column names, make a dict with the ascii names as keys instead #
    content = dict((corr[x], row[x]) for x in row.index if x in corr and row[x] is not numpy.nan)

    # Skip the case where the sample is not used #
    if content.get('used') == 'no': continue

    # Double quotes are used in JSON, replace them, and add quotes around #
    for k,v in content.items(): content[k] = '"' + unicode(v).replace('"', "'") + '"'

    # Special case for missing second contact #
    second_contact = {"contact_two_function": content['contact_two_function'],
                      "contact_two_name":     content['contact_two_name'],
                      "contact_two_email":    content['contact_two_email'],
                     } if content.get('contact_two_name') else False
    content['second_contact'] = second_contact

    # Figure out the path #
    path = home + "repos/sifes/metadata/json/projects/%s/%s/%s.json"
    path = path % (content['organization'].strip('"'),
                   content['project_short_name'].strip('"'),
                   content['sample_short_name'].strip('"'))

    # Create directory if it doesn't exist #
    dir_path = os.path.dirname(path)
    if not os.path.exists(dir_path): os.makedirs(dir_path)

    # Get the template and a renderer #
    template = JsonTemplate(content)
    renderer = pystache.Renderer(escape          = lambda u: u,
                                 string_encoding = 'utf8',
                                 file_encoding   = 'utf8',)

    # Make the text from the template, the content and a renderer #
    text = renderer.render(template, content)

    # Write it #
    with codecs.open(path, 'w', encoding='utf-8') as handle: handle.write(text)
