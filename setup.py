from distutils.core import setup

setup(
      name             = 'gefes',
      version          = '0.0.2',
      description      = 'Genome Extraction From Environmental Sequencing',
      long_description = open('README.txt').read(),
      license          = 'MIT',
      url              = 'http://github.com/limno/gefes/',
      author           = 'Lucas Sinclair',
      author_email     = 'lucas.sinclair@me.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['gefes'],
      requires         = ['sh', 'decorator', 'biopython', 'threadpool', 'patsy', 'scipy', 'matplotlib', 'statsmodels', 'pandas', 'ipython', 'scikit-learn', 'fastqident', 'rpy2'],
)


