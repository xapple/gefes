from distutils.core import setup

setup(
      name             = 'gefes',
      version          = '0.0.1',
      description      = 'Genome Extraction From Environmental Sequencing',
      long_description = open('README.txt').read(),
      license          = 'MIT',
      url              = 'http://github.com/limno/gefes/',
      author           = 'Lucas Sinclair',
      author_email     = 'lucas.sinclair@me.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['gefes'],
      requires         = ['biopython', 'matplotlib', 'sh', 'pandas', 'statsmodels', 'threadpool'],
)