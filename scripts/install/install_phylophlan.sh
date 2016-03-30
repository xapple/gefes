# Install old version from 2013 #
cd ~/test/
wget https://bitbucket.org/nsegata/phylophlan/get/default.zip
unzip default.zip
rm default.zip
mv *-phylophlan-* phylophlan
cd phylophlan
echo "TODO"

# Install 1.1.0 version from 2014 #
cd ~/programs/
hg clone https://bitbucket.org/nsegata/phylophlan
cd phylophlan
hg update 1.1.0
sed -i "s/usearch/usearch5/g" $HOME/programs/phylophlan/phylophlan.py
escaped_home=${HOME//\//\\\/}
sed -i "s/taxcuration\/')/$escaped_home\/programs\/phylophlan\/taxcuration\/')/g" $HOME/programs/phylophlan/phylophlan.py
export PATH="$HOME/programs/phylophlan":$PATH