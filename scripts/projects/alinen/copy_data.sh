# Run004 raw from sisu #
rsync -e 'ssh -o PreferredAuthentications=keyboard-interactive -o PubkeyAuthentication=no' -avz --progress --delete --copy-unsafe-links --checksum bob@milou.uppmax.uu.se:/home/bob/proj83/INBOX/140612_D00457_0037_AC49JKACXX $HOME/GEFES/raw/140612_D00457_0037_AC49JKACXX

# Run004 raw from uppmax #
rsync -r --progress ~/proj83/INBOX/140612_D00457_0037_AC49JKACXX/ bob@sisu-login1.csc.fi:/homeappl/home/bob/GEFES/raw/140612_D00457_0037_AC49JKACXX/

# Run004 raw link #
ln -s $HOME/GEFES/raw/140612_D00457_0037_AC49JKACXX $HOME/proj/b2014083/INBOX/140612_D00457_0037_AC49JKACXX

# Alinen cleaned #
rsync -r --progress ~/GEFES/views/samples/run004-* bob@sisu-login1.csc.fi:~/GEFES/views/samples/
