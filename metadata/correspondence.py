#!/usr/bin/env python2
# -*- coding: utf-8 -*-

###############################################################################
corr = {
    u'Project short name\n(no spaces and only ascii)': 'project_short_name',
    u'Project long name (free text)':                  'project_long_name',
    u'Sample short name\n(no spaces and only ascii)':  'sample_short_name',
    u'Sample long name (free text)':                   'sample_long_name',
    u'Sample #':                                       'sample_num',
    u'Used':                                           'used',

    u'Organization':                         'organization',
    u'Prefix':                               'prefix',
    u'Directory':                            'directory',
    u'Suffix':                               'suffix',

    u'Contact 1 function':                   'contact_one_function',
    u'Contact 1 name':                       'contact_one_name',
    u'Contact 1 email':                      'contact_one_email',
    u'Contact 2 function':                   'contact_two_function',
    u'Contact 2 name':                       'contact_two_name',
    u'Contact 2 email':                      'contact_two_email',

    u'Custom barcode':                       'custom_barcode',
    u'Multiplexed in':                       'multiplexed_in',
    u'Pool barcode':                         'pool_barcode',

    u'Header barcode':                       'header_barcode',
    u'Header rev. comp.':                    'header_rev_comp',
    u'Multiplex group reference':            'multiplex_group',

    u'Forward index #':                      'forward_num',
    u'Forward index sequence':               'forward_mid',
    u'Reverse index #':                      'reverse_num',
    u'Reverse index sequence':               'reverse_mid',

    u'Barcode ref.':                         'barcode_num',
    u'DNA con. [ng/µl]':                     'dna_after',
    u'PhiX spiking':                         'phix_spiking',
    u'Forward filename':                     'fwd_filename',
    u'Reverse filename':                     'rev_filename',
    u'Forward reads count':                  'fwd_read_count',
    u'Reverse reads count':                  'rev_read_count',
    u'Forward MD5 checksum':                 'fwd_md5',
    u'Reverse MD5 checksum':                 'rev_md5',

    u'Library strategy':                     'library_strategy',
    u'Library source':                       'library_source',
    u'Library selection':                    'library_selection',
    u'Library layout':                       'library_layout',
    u'Platform':                             'platform',
    u'Instrument model':                     'instrument_model',
    u'Instrument software':                  'instrument_software',
    u'Forward read length':                  'forward_read_length',
    u'Reverse read length':                  'reverse_read_length',

    u'Organism':                             'organism',
    u'Environement biome (env_biome)':       'env_biome',
    u'Environement feature (env_feature)':   'env_feature',
    u'Environement material (env_material)': 'env_material',

    u'Sampling Date (YYYY-MM-DD)':           'date',
    u'Latitude (N)':                         'latitude',
    u'Longitude (E)':                        'longitude',
    u'Country':                              'country',
    u'Location (free text)':                 'location',
    u'Design description (free text)':       'design_description',

    u'Bioproject':                           'bioproject',
    u'Biosample':                            'biosample',
}

###############################################################################
extras = {
    u'Depth [m]':                            'lorem',
    u'pH':                                   'lorem',
    u'Temperature [℃]':                      'lorem',
    u'Oxygen (O2) [mg/L]':                   'lorem',
    u'Carbon dioxide (CO2) [µM]':            'lorem',
    u'Methane (CH4) [µM]':                   'lorem',
    u'Iron II [µM]':                         'lorem',
    u'Iron III [µM]':                        'lorem',
    u'TOC [mg/L]':                           'lorem',
    u'SUVA [mg/L*m]':                        'lorem',
    u'total P [µg/L]':                       'lorem',
    u'total N [µg/L]':                       'lorem',
    u'Conductivity [μS/cm]':                 'lorem',
    u'Cell counts [cells/mL]':               'lorem',
    u'Secchi depth [cm]':                    'lorem',
    u'TOP [µg/L]':                           'lorem',
}