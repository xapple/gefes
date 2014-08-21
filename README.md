## GEFES

The acronym "GEFES" stands for **G**​enome **E**​xtraction **F**​rom **E**​nvironmental **S**​equencing.

This is yet another pipeline for assembling the short reads produced by shotgun-metagenomic sequencing experiments in an attempt to recompose full microbial genomes. With this tool, we could like to reconstitute the functional potential of the important bacterial and archaeal players in aquatic environments.

### Warnings:

* First of all, this is still very much work in progress. We are not yet at a stage where full metabolic predictions can be made.

* Secondly, of course these kind of approaches are not entirely revolutionary. Other labs on the globe have been doing things like this for a two or three years. But the main reason we are replicating these efforts in a way is to acquire our own setup that we are able to control. Other labs do make their tools available, but if you manage to install them on your computer you end up with something that was tailored for their type of data and that you can't change. Moreover, you you don't really understand how it works. With our own setup we have a flexibility that let's us try new ideas and add them in very quickly.

* You also have to keep in mind that what we are measuring is a potential. And it's not because you see a particular gene in a population or clade that it is always being actively translated, transcribed, folded and exported. Ideally, you should use the results from such tools like GEFES to direct your laboratory assays which will truly confirm if a given process is taking place.

### Context:

As you know almost all microbes can't be isolated or cultured easily. So instead, we go and do metagenomic sampling. We take a glass of water from a lake, do a total DNA extraction and insert the solution into a sequencer. As a result we get a file full of short DNA reads each coming from a different microbe.

It's quite different from when you are able to isolate a bacterium such as E. Coli. In that case your DNA reads are coming from random locations of the chromosome, but they are all from a copy of the same genome. This makes it easy to pieces things together afterwards.

What we have as starting data in our case is more messy. Every DNA read is potentially coming from a different species. Plus, the fragments are really short and only span a fraction of a typical microbial gene.

How do we put the reads together to make genomes ? How are we going to figure out which short sequence was coming from which species ? This is what GEFES is supposed to help with.

### Overview:

Starting from the raw reads there are about eight distinct processing steps in GEFES:

1. Cleaning
2. Co-Assemble
3. Map
4. Cluster contigs
5. Validate
6. Assign taxonomy
7. Annotate
8. Predict metabolism

Unfortunately, no other detailed documentation has been written yet but the code is clean and commented. In addition these two descriptive files might help you figure out what is going on:

* documentation/diagram.pdf
* documentation/flowchart.pdf

## Installing

No automated installation has been developed for the GEFES package.
But following this document and typing these commands on your bash prompt should get you started. It is designed so you don't need super user privileges at any step.
If you cannot get a functional installation set up, contact the authors.

### Step 1: Cloning the repository
Here you will download a copy of the code from github and place it in your home directory.

    $ cd ~
    $ mkdir repos
    $ cd repos
    $ git clone git@github.com:limno/gefes.git

### Step 2A: Modify your python search path
Here you will edit your ``.bashrc`` or ``.bash_profile`` to add a reference to the module you just downloaded. When you type `import gefes` python will know where to look !

    $ vim ~/.bash_profile
    export PYTHONPATH="$HOME/repos/gefes/":$PYTHONPATH

### Step 2B: Modify your binary search path
There are some dependencies that are not too large and can be bundled as binaries with this module. To avoid the hassle of downloading and compiling these requirements, we have added a bunch of them in the repository. Add them to your `$PATH`.

    $ vim ~/.bash_profile
    export PATH="$HOME/repos/gefes/bin/linux/:$PATH"
    export PATH="$HOME/repos/gefes/bin/linux/mcl/:$PATH"

### Step 3: Install your own version of python
Your system probably comes with a version of python installed. But the variations from system to system are too great to rely on any available python. We prefer to just install our own in the home directory.

For this we will be using this excellent project: https://github.com/yyuu/pyenv

To install it you may use this sister project: https://github.com/yyuu/pyenv-installer

Basically you just need to type this command:

    $ curl -L https://raw.githubusercontent.com/yyuu/pyenv-installer/master/bin/pyenv-installer | bash

These lines go into your ``.bash_profile``:

    $ vim ~/.bash_profile
    export PYENV_ROOT="$HOME/.pyenv"
    export PATH="$PYENV_ROOT/bin:$PATH"
    eval "$(pyenv init -)"

Relaunch your shell and type these commands to get the right version of python now:

    pyenv install 2.7.6
    pyenv rehash
    pyenv global 2.7.6

### Step 4: Install all required python packages
GEFES uses many third party python libraries. You can get them by running these commands:

    $ pip install sh
    $ pip install decorator
    $ pip install biopython
    $ pip install threadpool
    $ pip install patsy
    $ pip install scipy
    $ pip install matplotlib
    $ pip install pandas
    $ pip install statsmodels
    $ pip install ipython
    $ pip install scikit-learn
    $ pip install fastqident
    $ pip install rpy2
    $ pip install pysam

Don't forget to rehash the binary links at the end:

    $ pyenv rehash

### Step 5: Make a working directory with the raw data linked
GEFES will search for the raw reads in a directory called ``INBOX`` somewhere in your home. You should set up the correct symbolic links (or modify the code so that it searches for the FASTQ files where you want it to).

    $ cd ~
    $ mkdir GEFES
    $ cd GEFES
    $ mkdir views
    $ mkdir INBOX
    $ cd INBOX
    $ ln -s /proj/b2011035/INBOX/* ./

### Step 6: Obtaining extra dependencies
GEFES also makes use of many third party programs which need to be installed and accessible from your ``$PATH``. These dependencies each have specific installation procedures and include:

 * [sickle](https://github.com/najoshi/sickle) version X.X.X providing ``sickle``
 * [Ray](http://sourceforge.net/projects/denovoassembler/) version X.X.X providing ``Ray23``
 * [htslib](https://github.com/samtools/htslib) (needed by samtools) version X.X.X
 * [samtools](http://samtools.sourceforge.net) version X.X.X providing ``samtools``
 * [bedtools]() version X.X.X providing ``genomeCoverageBed``
 * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) version X.X.X providing ``bowtie2``
 * [Glimmer]() version X.X.X providing ``XXX``
 * [Prodigal]() version X.X.X providing ``XXX``
 * [Spades]() version X.X.X providing ``XXX``
 * get-motif-counts.awk
 * long-orfs
 * extract
 * build-icm

These can take some time to install and unfortunately we can't package them with our project !

### Step 7: Troubleshooting
You might run into difficulties at some of the steps above -- if that is the case check this section for a solution.

#### - "pip install scipy" missing BLAS:
When you install scipy, you might need these two dependencies before hand if they are not on your system: http://stackoverflow.com/questions/7496547/python-scipy-needs-blas

#### - git version too old:
If your git version is too old you can follow these steps: http://blog.justin.kelly.org.au/install-git-in-your-home-directory-how-to/

#### - Compiling Ray:
To compile Ray on a compute cluster you might have to do something this:

    $ module swap PrgEnv-intel PrgEnv-gnu
    $ make -j8 MPI_IO=y MPICXX=cc MAXKMERLENGTH=91

Or something like this:

    $ module load openmpi/1.4.5
    $ module load pgi
    $ make -j8 MPI_IO=y MPICXX=mpicc MAXKMERLENGTH=91
