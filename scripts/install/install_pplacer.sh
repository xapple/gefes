cd ~/test/
wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha17/pplacer-Linux-v1.1.alpha17.zip
unzip pplacer-Linux-v1.1.alpha17.zip
mv pplacer-Linux-v1.1.alpha17 $HOME/programs/pplacer17/
export PATH="$HOME/programs/pplacer17":$PATH
find $HOME/programs/pplacer17 -type f -print0 | xargs -0 chmod u-w