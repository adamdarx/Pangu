# wget https://github.com/doreamon-design/clash/releases/download/v2.0.24/clash_2.0.24_linux_amd64.tar.gz
# tar zxvf clash_2.0.24_linux_amd64.tar.gz
# chmod +x clash
# mv clash /usr/local/bin/clash
# mkdir /etc/clash
# cat > /etc/clash/config.yaml << EOF
# # Clash 配置示例
# # 直接从客户端软件复制
# EOF
# echo "export https_proxy=http://127.0.0.1:7890 http_proxy=http://127.0.0.1:7890 all_proxy=socks5://127.0.0.1:7891" >> ~/.bashrc
# source ~/.bashrc

wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.10.tar.gz
tar -zxvf openmpi-5.0.10.tar.gz
cd openmpi-5.0.10
./configure --prefix=/openmpi --with-cuda=/usr/lib/cuda --enable-mpi-cxx --enable-shared
make all install
export PATH=$PATH:/openmpi/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/openmpi/lib
source ~/.bashrc

wget https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_3/downloads/hdf5-1.14.3.tar.gz
tar -zxvf hdf5-1.14.3.tar.gz
export CC=mpicc
export HDF5_MPI="ON"
./configure --enable-shared --enable-parallel --prefix=/hdf5
export HDF5_DIR="/hdf5"
make
make install
make check-install
export PATH=$PATH:"/hdf5/bin"
export HDF5_ROOT="/hdf5"
source ~/.bashrc