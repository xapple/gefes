cd ~/test/
wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
tar -xjf samtools-1.2.tar.bz2
rm samtools-1.2.tar.bz2
cd samtools-1.2

vim Makefile
# LIBCURSES=  -lcurses
# change to
# LIBCURSES=  -lncurses

mkdir $HOME/programs/samtools/
make prefix=$HOME/programs/samtools install
export PATH="$HOME/programs/samtools/bin":$PATH