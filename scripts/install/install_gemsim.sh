cd ~/test/
wget http://downloads.sourceforge.net/project/gemsim/GemSIM_v1.6.tar.gz
tar -xzf GemSIM_v1.6.tar.gz

chmod +x GemSIM_v1.6/*.py
mv GemSIM_v1.6 ../programs/gemsim

export PATH="$HOME/programs/gemsim":$PATH
