# Prokka #
cd ~/test/
wget http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift/phylosift_latest.tar.bz2
tar xjf phylosift_latest.tar.bz2
rm phylosift_latest.tar.bz2
mv phylosift_* ~/programs/phylosift
export PATH="$HOME/programs/phylosift":$PATH
