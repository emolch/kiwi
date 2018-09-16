#!/bin/bash

set -e

hdf5=hdf5-1.8.21

if [ ! -f "${hdf5}.tar.gz" ] ; then
    wget "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.21/src/${hdf5}.tar.gz"
fi

if [ -d "${hdf5}" ] ; then 
    rm -rf "${hdf5}"
fi

tar -xzf "${hdf5}.tar.gz"

cd "$hdf5"
curdir=`pwd`
./configure --enable-fortran --prefix="$curdir/hdf5"
make -j 2
make install

echo "set these variables in Makefile.local:"
echo
echo "  HDF5DIR = hdf5-1.8.21/hdf5"
echo "  INCHDF = -I$(HDF5DIR)/include"
echo "  LIBHDF = $(HDF5DIR)/lib/libhdf5_fortran.a $(HDF5DIR)/lib/libhdf5.a -ldl -lz"
echo
