## github actions版本
export NETCDF="/opt/netcdf_v4.9.2_openmpi"
export PNDIR="$NETCDF"
export LD_LIBRARY_PATH="$NETCDF/lib"

# export PATH=$NETCDF/bin:$PATH
# export NETCDF=$NETCDF
export NETCDF_C_PATH=$NETCDF
export NETCDF_FORTRAN_PATH=$NETCDF
export PNETCDF_PATH=$PNDIR

export CFLAGS="-fPIC -shared"
export FFLAGS="-fPIC -shared"

# User-specific variables
BUILD_DIR="bld"
INSTALL_DIR="eclm"

mkdir bld

# Run cmake
cmake -S src -B "$BUILD_DIR" \
      -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz -fPIC" \
      -DCMAKE_Fortran_COMPILER=mpifort

cd bld
make -j8
