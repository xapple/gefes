## gefes version 0.0.5

The acronym `gefes` stands for **G**​enome **E**​xtraction **F**​rom **E**​nvironmental **S**​equencing.

This is yet another pipeline for assembling the short reads produced by shotgun-metagenomic sequencing experiments in an attempt to recompose full microbial genomes. With this tool, we could like to reconstitute the functional potential of the important bacterial and archaeal players in aquatic environments.

It was developed by Lucas Sinclair (<lucas.sinclair@me.com>) while working in the Limnology department at the Evolution Biology Center of Uppsala University. The code has an MIT license and everyone is welcome to use, modify or extend the pipeline.

## Warnings:

* First of all, this is still very much work in progress. We are not yet at a stage where full metabolic predictions can be made.

* Secondly, of course these kind of approaches are not entirely revolutionary. Other labs on the globe have been doing things like this for a two or three years. But the main reason we are replicating these efforts in a way is to acquire our own setup that we are able to control. Other labs do make their tools available, but if you manage to install them on your computer you end up with something that was tailored for their type of data and that you can't change. Moreover, you you don't really understand how it works. With our own setup we have a flexibility that lets us try new ideas and implement them in the code very quickly.

* Thirdly, you also have to keep in mind that what we are measuring is a potential. It's not because you see a particular gene in a population or strain that it is necessarily being actively translated, transcribed, folded and exported. Ideally, you should use the results from such tools like `gefes` to direct your laboratory assays which will truly confirm if a given process is taking place.

* Finally, the `gefes` project is not a biologist-oriented tool that supports all the possible use cases one could have with metagenomic sequence data out of the box. For instance, it does not have a graphical interface to operate, nor any bash/sh/csh commands. Indeed, as each sequencing experiment will have different goals and scientific questions associated to it, there cannot be a standard set of procedures to apply to the data. Instead, the `gefes` project a flexible and modular collections of packages written in proper, clean and commented object-oriented python which enables the user to survey, modify and extend the code-base easily -- provided he has a sufficient knowledge in programming. It is a basis upon which the scientist can set up the processing and analysis that he sees fit for his own data sparing him from having to develop lots of the infrastructure needed himself.

## Context:

As you know almost all microbes living in natural environments can't be isolated or cultured easily. So instead, we go and do shotgun metagenomic sampling. We take a glass of water from a lake, do a total DNA extraction and insert the solution into a high-throughput sequencer. As a result, we receive a file full of short DNA reads each coming (statistically) from a different microbe.

It's quite different from when you are able to isolate a bacterium such as E. Coli. In that case, your DNA reads are coming from random locations of the chromosome, but they are all originating from a copy of the same genome. This makes it easy to pieces things together afterwards.

What we have as starting data in our case is more messy. Every DNA read is potentially coming from a different species. Plus, the fragments are really short and only span a fraction of a typical microbial gene.

How do we put the reads together to make genomes ? How are we going to figure out which short sequence was coming from which species ? This is what `gefes` is supposed to help with.

Many objects common to any analysis such as a "FASTQ file pair", a "Sample", a "Aggregate of Samples", a "Sequence quality checker", an "Assembly", a "Read mapper" are provided. In addition you will find routines for sending these objects through well-known algorithms such as Sickle, Ray, Bowtie, etc. Lots of other functionality is also present such as a multitude of visualization in `matplotlib` and other things such as the ability to automatically distribute the computation on a network of computers (via SLURM).

## Overview:

Starting from the raw reads there are about eight distinct processing steps in `gefes`:

1. Cleaning
2. Co-Assemble
3. Map
4. Cluster contigs
5. Validate
6. Assign taxonomy
7. Annotate
8. Predict metabolism

Unfortunately, no other detailed documentation has been written yet but the code is clean and commented. In addition, these two descriptive files might help you figure out what is going on:

* [Flowchart of data processing](/../master/documentation/flowchart.pdf?raw=true "Flowchart")
* [Tentative UML diagram of objects compositionn](/../master/documentation/diagram.pdf?raw=true "Diagram")

## Installing

No automated installation has been developed for the `gefes` package. Ultimately you will be able to install it like this:

~~~
To install `seqenv` onto your machine, use the python package manager:

    $ pip install gefes

You might be installing this onto a computer server which you don't own and thus don't have sufficient privileges. In that case you can install everything in your home directory like this:

    $ pip install --user gefes

If this still doesn't work, you might be missing the `pip` program on your system or the correct version of Python (any version `2.7.x`). You can get both of these things by using using this little project: https://github.com/yyuu/pyenv
~~~

In the meantime, following this document and typing these commands on your bash prompt should get you started. It is designed so you don't need super user privileges at any step. If you cannot get a functional installation set up, contact the authors.

#### Step 1: Cloning the repository
Here you will download a copy of the code from github and place it in your home directory.

    $ cd ~
    $ mkdir repos
    $ cd repos
    $ git clone git@github.com:limno/gefes.git

#### Step 2A: Modify your python search path
Here you will edit your ``.bashrc`` or ``.bash_profile`` to add a reference to the module you just downloaded. When you type `import gefes` python will know where to look !

    $ vim ~/.bash_profile
    export PYTHONPATH="$HOME/repos/gefes/":$PYTHONPATH

#### Step 2B: Modify your binary search path
There are some dependencies that are not too large and can be bundled as binaries with this module. To avoid the hassle of downloading and compiling these requirements, we have added a bunch of them in the repository. Add them to your `$PATH`.

    $ vim ~/.bash_profile
    export PATH="$HOME/repos/gefes/bin/linux/:$PATH"
    export PATH="$HOME/repos/gefes/bin/linux/mcl/:$PATH"

#### Step 3: Install your own version of python
Your system probably comes with a version of python installed. But the variations from system to system are too great to rely on any available setup. We strongly recommand to just install our own in the home directory. Also, we will then be able to install modules without administrator privileges. You can skip this step if you are confident enough about what you are doing.

For this we will be using this excellent project: https://github.com/yyuu/pyenv

To install it you may use this sister project: https://github.com/yyuu/pyenv-installer

Basically you just need to type this command:

    $ curl -L https://raw.githubusercontent.com/yyuu/pyenv-installer/master/bin/pyenv-installer | bash

These lines go into your ``.bash_profile``:

    $ vim ~/.bash_profile
    export PYENV_ROOT="$HOME/.pyenv"
    export PATH="$PYENV_ROOT/bin:$PATH"
    eval "$(pyenv init -)"

Relaunch your shell and type these commands to get the right version of python:

    $ pyenv install 2.7.8
    $ pyenv rehash
    $ pyenv global 2.7.8

#### Step 4: Install all required python packages
`gefes` uses many third party python libraries. You can get them by running these commands:

    $ pip install plumbing
    $ pip install fasta
    $ pip install pymarktex
    $ pip install matplotlib
    $ pip install pandas
    $ pip install biopython
    $ pip install sh
    $ pip install tqdm
    $ pip install decorator
    $ pip install shell_command
    $ pip install threadpool
    $ pip install patsy
    $ pip install scipy
    $ pip install statsmodels
    $ pip install scikit-learn
    $ pip install pysam

Don't forget to rehash the executable links at the end if you are using pyenv like we do:

    $ pyenv rehash

#### Step 5: Obtaining extra dependencies
`gefes` also makes use of many third party programs which need to be installed and accessible from your ``$PATH``. Depending on what parts of the pipeline you are planning to run, you might not need them all. You can try and install the missing ones when `gefes` complains about a missing executable. These dependencies each have specific installation procedures and include:

 * [Ray](http://sourceforge.net/projects/denovoassembler/) version 2.3.1 providing ``ray231`` (included in repository)
 * [Picard Tools](http://broadinstitute.github.io/picard/) version 1.101 providing ``MarkDuplicates.jar`` (included in repository)
 * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) version 2.2.4 providing ``bowtie2`` and ``bowtie2_build2``
 * [bedtools](http://bedtools.readthedocs.org/en/latest/) version 2.15.0 providing ``genomeCoverageBed``
 * [samtools](http://samtools.sourceforge.net) version 0.1.19 providing ``samtools``

 * [Glimmer]() version X.X.X providing ``XXX``
 * [Prodigal]() version X.X.X providing ``XXX``
 * [Spades]() version X.X.X providing ``XXX``
 * long-orfs
 * extract
 * build-icm

These can take some time to install and unfortunately we can't package them with our project !

#### Step 6: Make a working directory with the raw data linked
`gefes` will generate results in a directory named `views` in a directory named `GEFES` in your home folder. You can change this by editing the code of course.

    $ cd ~
    $ mkdir GEFES
    $ cd GEFES
    $ mkdir views
    $ mkdir INBOX

`gefes` will usually search for the raw reads in a directory called `INBOX` somewhere in your home. You should set up the correct symbolic links.

    $ cd ~
    $ cd GEFES
    $ mkdir INBOX
    $ cd INBOX
    $ ln -s /proj/b2011035/INBOX/* ./

#### Step 7: Make a JSON file describing your project
`gefes` will create `Project` objects with associated `Sample` objects based on a user inputted JSON files. Look at the `json` directory at the root of the repository and make your own.

#### Step 8: Troubleshooting
You might run into difficulties at some of the steps above -- if that is the case check this section for a solution.

##### - "pip install scipy" missing BLAS:
When you install scipy, you might need these two dependencies before hand if they are not on your system: http://stackoverflow.com/questions/7496547/python-scipy-needs-blas

##### - git version too old:
If your git version is too old you can follow these steps: http://blog.justin.kelly.org.au/install-git-in-your-home-directory-how-to/

##### - Compiling Ray:
To compile Ray on a compute cluster you might have to do something this:

    $ module swap PrgEnv-intel PrgEnv-gnu
    $ make -j24 MPI_IO=y MPICXX=CC MAXKMERLENGTH=91

Or something like this:

    $ module load openmpi/1.4.5
    $ module load pgi
    $ make -j16 MPI_IO=y MPICXX=mpicc MAXKMERLENGTH=91

##### - "pip install matplotlib" missing freetype:
If you are on OS X you can simply fix this error by typing `$ brew install freetype`

## Acknowledgments
A special thanks to all those who helped create this pipeline and make it as great as it is:
* Ino De Bruijn [@inodb](https://github.com/inodb)
* Moritz Buck [@moritzbuck](https://github.com/moritzbuck)
* Umer Zeeshan Ijaz [@umerijaz](https://github.com/umerijaz)

## Flowchart
Below is drawn the flowchart describing the data processing along all the steps of GEFES:

![Flowchart](/../master/documentation/flowchart.png?raw=true "Flowchart")
