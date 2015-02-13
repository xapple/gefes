# Install #
cd ~/test/
wget http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift/phylosift_latest.tar.bz2
tar xjf phylosift_latest.tar.bz2
rm phylosift_latest.tar.bz2
mv phylosift_* ~/programs/phylosift
export PATH="$HOME/programs/phylosift":$PATH

# Databases #
mkdir -p ~/share/phylosift
cd ~/share/phylosift
wget http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift_markers/markers.tgz
wget http://edhar.genomecenter.ucdavis.edu/~koadman/ncbi.tgz
tar xzf markers.tgz
tar xzf ncbi.tgz
rm ncbi.tgz
rm markers.tgz

# Run a test example #
cd ~/repos/examples/phylosift
rm -rf results
phylosift all --output=results contig.fasta
