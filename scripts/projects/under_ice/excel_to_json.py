#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script that reads an excel file and produces our correctly formated JSON files.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ excel_to_json.py ~/repos/gefes/scripts/projects/under_ice/meta_data.xlsx
"""

# Modules #
import sys, os, pandas, codecs, numpy, datetime

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
template_proj = u"""{
    "contacts": {
        "%(contact_one_function)s": {
            "name": "%(contact_one_name)s",
            "email": "%(contact_one_email)s"
        },
        "%(contact_two_function)s": {
            "name": "%(contact_two_name)s",
            "email": "%(contact_two_email)s"
        }
    },

    "project":       "%(project_short_name)s",
    "project_name":  "%(project_long_name)s",
    "project_num":   %(project_num)s,

    "uppmax_project_id":    "%(uppmax_project_id)s",
    "illumina_run_id":      "%(illumina_run_id)s",

    "library_strategy":     "%(library_strategy)s",
    "library_source":       "%(library_source)s",
    "library_selection":    "%(library_selection)s",
    "library_layout":       "%(library_layout)s",
    "platform":             "%(platform)s",
    "instrument_model":     "%(instrument_model)s",
    "instrument_software":  "%(instrument_software)s",
    "forward_read_length":  %(forward_read_length)s,
    "reverse_read_length":  %(reverse_read_length)s,

    "date":                 "%(date)s",
    "latitude":             [%(latitude)s, "N"],
    "longitude":            [%(longitude)s, "E"],
    "country":              "%(country)s",
    "location":             "%(location)s",

    "bioproject":           "%(bioproject)s",

    "samples_base_dir":     "%(samples_base_dir)s",
    "samples": [

%(samples_json)s
    ]
}"""

###############################################################################
template_sample = u"""    {
        "sample_name":          "%(sample_short_name)s",
        "sample_long_name":     "%(sample_long_name)s",
        "sample_directory":     "%(sample_directory)s",
        "sample_num":           %(sample_num)s,

        "forward_reads":        "%(fwd_filename)s",
        "reverse_reads":        "%(rev_filename)s",
        "forward_md5":          "%(fwd_md5)s",
        "reverse_md5":          "%(rev_md5)s",
        "forward_read_count":   %(fwd_count)s,
        "reverse_read_count":   %(rev_count)s,

        "organism":             "%(organism)s",
        "env_biome":            "%(env_biome)s",
        "env_feature":          "%(env_feature)s",
        "env_material":         "%(env_material)s",
        "design_description":   "%(design_description)s",

        "biosample":            "%(biosample)s",

        "depth":                [%(depth)s, "m"],
        "ph":                   %(ph)s,
        "toc":                  [%(toc)s, "mg/l"],
        "ton":                  [%(ton)s, "mg/l"],
        "top":                  [%(top)s, "µg/l"],
        "sulfate":              [%(sulfate)s, "mg/l"],
        "oxygen":               [%(oxygen)s, "mg/l"],
        "conductance":          [%(conductance)s, "µS/cm"],
        "temperature":          [%(temperature)s, "Celsius"],
        "filtered_volume":      [%(filtered_volume)s, "ml"],
        "cell_counts":          [%(cell_counts)s, "cells/ml"],
        "co2":                  [%(co2)s, "µM"],
        "ch4":                  [%(ch4)s, "µM"],
        "feII":                 [%(fe2)s, "µM"],
        "feIII":                [%(fe3)s, "µM"],
        "fe_total":             [%(fe_total)s, "µM"],
        "suva":                 [%(suva)s, "mg/l*m"]
    }"""

###############################################################################
correspondence = {
    u'Uppmax ref':                           'uppmax_project_id',
    u'Run name':                             'illumina_run_id',
    u'Sample name':                          'sample_directory',
    u'Base directory':                       'samples_base_dir',

    u'Forward filename':                     'fwd_filename',
    u'Reverse filename':                     'rev_filename',
    u'Forward reads count':                  'fwd_count',
    u'Reverse reads count':                  'rev_count',
    u'Forward MD5 checksum':                 'fwd_md5',
    u'Reverse MD5 checksum':                 'rev_md5',

    u'Contact 1 function':                   'contact_one_function',
    u'Contact 1 name':                       'contact_one_name',
    u'Contact 1 email':                      'contact_one_email',
    u'Contact 2 function':                   'contact_two_function',
    u'Contact 2 name':                       'contact_two_name',
    u'Contact 2 email':                      'contact_two_email',

    u'Project short name (no spaces and only ascii)': 'project_short_name',
    u'Project long name (free text)':                 'project_long_name',
    u'Project #':                                     'project_num',
    u'Sample short name (no spaces and only ascii)':  'sample_short_name',
    u'Sample long name (free text)':                  'sample_long_name',
    u'Sample #':                                      'sample_num',

    u'Library strategy':             'library_strategy',
    u'Library source':               'library_source',
    u'Library selection':            'library_selection',
    u'Library layout':               'library_layout',
    u'Platform':                     'platform',
    u'Instrument model':             'instrument_model',
    u'Instrument software':          'instrument_software',
    u'Forward read length':          'forward_read_length',
    u'Reverse read length':          'reverse_read_length',

    u'Organism':                     'organism',
    u'Environement biome':           'env_biome',
    u'Environement feature':         'env_feature',
    u'Environement material':        'env_material',

    u'Sampling Date (YYYY-MM-DD)':   'date',
    u'Latitude (N)':                 'latitude',
    u'Longitude (E)':                'longitude',
    u'Country':                      'country',
    u'Location (free text)':         'location',
    u'Design description':           'design_description',

    u'Bioproject':                   'bioproject',
    u'Biosample':                    'biosample',

    u'Depth':                        "depth",
    u'pH':                           "ph",
    u"TOC average":                  "toc",
    u"TON average":                  "ton",
    u"TOP average":                  "top",
    u"Sulfate average":              "sulfate",
    u"Oxygen":                       "oxygen",
    u"Cond.":                        "conductance",
    u"Temp.":                        "temperature",
    u"Volume chosen":                "filtered_volume",
    u"Cell counts":                  "cell_counts",
    u"CO2":                          "co2",
    u"CH4":                          "ch4",
    u"FeII":                         "fe2",
    u"FeIII":                        "fe3",
    u"Fe Total":                     "fe_total",
    u"SUVA":                         "suva",
}

###############################################################################
# Get the shell arguments #
if len(sys.argv) < 2: sys.exit(sys.modules[__name__].__doc__)
xlsx_path = sys.argv[1]

# Check that the path is valid #
if not os.path.exists(xlsx_path): raise Exception("No file at %s." % xlsx_path)

# Load data #
xlsx = pandas.ExcelFile(xlsx_path)
df = xlsx.parse(xlsx.sheet_names[0])

# Identify projects #
proj_header = 'Project short name (no spaces and only ascii)'
projects = set(df[proj_header]) - set((numpy.nan,))

# Make as many JSON files as projects #
for proj in projects:
    # Identify rows #
    rows = df.loc[df[proj_header] == proj]
    # Now make each sample #
    samples_text = []
    for i, row in rows.iterrows():
        sample_data = dict((correspondence[x], row[x]) for x in row.index if x in correspondence)
        # Skip samples missing short names #
        if sample_data.get('sample_short_name') is numpy.nan: continue
        # Sometimes the time is introduced in the date #
        if isinstance(sample_data['date'], datetime.datetime): sample_data['date'] = sample_data['date'].date().isoformat()
        # Sometimes the sample number becomes a float #
        sample_data['sample_num'] = int(sample_data['sample_num'])
        # Template the text #
        text = template_sample % sample_data
        text = text.replace('"nan",', 'null,')
        text = text.replace('nan,',   'null,')
        samples_text.append(text)
    samples_text = ',\n\n'.join(samples_text)
    # Use the first row for project meta data #
    row = rows.iloc[0]
    project_data = dict((correspondence[x], row[x]) for x in row.index if x in correspondence)
    project_data['samples_json'] = samples_text
    # Sometimes the time is introduced in the date #
    if isinstance(project_data['date'], datetime.datetime): project_data['date'] = project_data['date'].date().isoformat()
    # Sometimes the sample number becomes a float #
    project_data['project_num'] = int(project_data['project_num'])
    # Template the text #
    text = template_proj % project_data
    text = text.replace('"nan",', 'null,')
    text = text.replace('nan,',   'null,')
    # Figure out the path #
    path = home + "repos/gefes/json/%s.json" % project_data['project_short_name']
    # Write #
    with codecs.open(path, 'w', encoding='utf-8') as handle: handle.write(text)