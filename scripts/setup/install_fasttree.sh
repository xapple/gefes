cd ~/test/
wget http://www.microbesonline.org/fasttree/FastTree
mkdir $HOME/programs/fast_tree/
chmod +x FastTree
mv FastTree $HOME/programs/fast_tree/
export PATH="$HOME/programs/fast_tree":$PATH
