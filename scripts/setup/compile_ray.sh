cd ~/test
wget https://github.com/sebhtml/ray/archive/v2.3.1.zip
unzip v2.3.1
wget https://github.com/sebhtml/RayPlatform/archive/v2.0.1.zip
unzip v2.0.1
mv RayPlatform-2.0.1 RayPlatform
cd RayPlatform
module swap PrgEnv-cray PrgEnv-intel
make -j24 MPI_IO=y MPICXX=CC MAXKMERLENGTH=91
cp Ray ~/repos/gefes/bin/cray/ray231