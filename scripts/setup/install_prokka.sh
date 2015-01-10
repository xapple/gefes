# Dependencies #
mkdir -p ~/programs/cpanm/
cd ~/programs/cpanm/
curl -LO http://xrl.us/cpanm
chmod +x cpanm
export PATH="$HOME/programs/cpanm":$PATH
cpanm install --force Bio::Perl

# Prokka #
cd ~/test/
wget https://github.com/Victorian-Bioinformatics-Consortium/prokka/archive/v1.10.tar.gz
tar xzf v1.10.tar.gz
rm v1.10.tar.gz
mv prokka-1.10 ~/programs/prokka110/
cd ~/programs/prokka110/
export PATH="$HOME/programs/prokka110/bin":$PATH
prokka --setupdb