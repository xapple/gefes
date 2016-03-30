cd ~/test/
wget https://github.com/hyattpd/Prodigal/archive/v2.6.2.tar.gz
tar xzf v2.6.2.tar.gz
cd Prodigal-2.6.2
mkdir -p ~/programs/prodigal262/
make install INSTALLDIR=~/programs/prodigal262/
export PATH="$HOME/programs/prodigal262":$PATH
