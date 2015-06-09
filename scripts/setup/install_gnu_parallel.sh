cd ~/test
wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
tar -xjf parallel-latest.tar.bz2
cd parallel-*
./configure --prefix=$HOME && make && make install