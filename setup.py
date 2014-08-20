from distutils.core import setup

setup(
      name             = 'gefes',
      version          = '0.0.4',
      description      = 'Genome Extraction From Environmental Sequencing',
      long_description = open('README.md').read(),
      license          = 'MIT',
      url              = 'http://github.com/limno/gefes/',
      author           = 'Lucas Sinclair',
      author_email     = 'lucas.sinclair@me.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['gefes'],
      requires         = ['plumbing', 'fasta', 'sh', 'biopython', 'decorator', 'threadpool', 'patsy', 'scipy', 'matplotlib', 'statsmodels', 'pandas', 'ipython', 'scikit-learn', 'fastqident', 'rpy2'],
)