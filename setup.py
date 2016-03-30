from distutils.core import setup

setup(
      name             = 'gefes',
      version          = '0.1.5',
      description      = 'Genome Extraction From Environmental Sequencing',
      long_description = open('README.md').read(),
      license          = 'Proprietary software',
      url              = 'http://github.com/xapple/gefes/',
      author           = 'Lucas Sinclair',
      author_email     = 'lucas.sinclair@me.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['gefes'],
      requires         = ['plumbing', 'fasta', 'pymarktex', 'seqsearch', 'sh', 'tqdm', 'biopython', 'decorator', 'threadpool', 'scipy', 'matplotlib', 'pandas', 'ipython'],
)