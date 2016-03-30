# Dependencies #
cd ~/test/
wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz
tar -xzf jellyfish-1.1.11.tar.gz
rm jellyfish-1.1.11.tar.gz
cd jellyfish-1.1.11
mkdir $HOME/programs/jellyfish1/
./configure --prefix=$HOME/programs/jellyfish1/
make
make install
export PATH="$HOME/programs/jellyfish1/bin:$PATH"

# Install #
cd ~/test/
wget http://ccb.jhu.edu/software/kraken/dl/kraken-0.10.5-beta.tgz
tar -xzf kraken-0.10.5-beta.tgz
rm kraken-0.10.5-beta.tgz
cd kraken-0.10.5-beta
mkdir $HOME/programs/kraken/
./install_kraken.sh $HOME/programs/kraken/

# Link executables #
mkdir $HOME/programs/kraken/bin
ln -s $HOME/programs/kraken/kraken            $HOME/programs/kraken/bin/kraken
ln -s $HOME/programs/kraken/kraken-build      $HOME/programs/kraken/bin/kraken-build
ln -s $HOME/programs/kraken/kraken-filter     $HOME/programs/kraken/bin/kraken-filter
ln -s $HOME/programs/kraken/kraken-mpa-report $HOME/programs/kraken/bin/kraken-mpa-report
ln -s $HOME/programs/kraken/kraken-report     $HOME/programs/kraken/bin/kraken-report
export PATH="$HOME/programs/kraken/bin:$PATH"

# Build DB #
kraken-build --standard --threads 16 --db $HOME/databases/kraken/standard
