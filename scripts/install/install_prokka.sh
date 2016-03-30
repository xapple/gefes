# Dependencies #
mkdir -p ~/programs/cpanm/
cd ~/programs/cpanm/
curl -LO http://xrl.us/cpanm
chmod +x cpanm
export PATH="$HOME/programs/cpanm":$PATH
cpanm install --force Bio::Perl

# Prokka 110 #
cd ~/test/
wget http://www.vicbioinformatics.com/prokka-1.10.tar.gz
tar xzf prokka-1.10.tar.gz
rm prokka-1.10.tar.gz
mv prokka-1.10 ~/programs/prokka110/
cd ~/programs/prokka110/
export PATH="$HOME/programs/prokka110/bin":$PATH
prokka --setupdb

# Prokka 111 #
cd ~/test/
wget http://www.vicbioinformatics.com/prokka-1.11.tar.gz
tar xzf prokka-1.11.tar.gz
rm prokka-1.11.tar.gz
mv prokka-1.11 ~/programs/prokka111/
export PATH="$HOME/programs/prokka111/bin":$PATH
prokka --setupdb

# Fix a problem #
wget https://github.com/tseemann/prokka/blob/master/binaries/linux/tbl2asn?raw=true
mv tbl2asn\?raw\=true tbl2asn
chomd +x tbl2asn
mv tbl2asn ~/programs/prokka111/binaries/linux/tbl2asn