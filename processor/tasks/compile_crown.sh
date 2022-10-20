#! /bin/bash
CROWNFOLDER=$1
ANALYSIS=$2
CONFIG=$3
SAMPLES=$4
ERAS=$5
SCOPE=$6
SHIFTS=$7
INSTALLDIR=$8
BUILDDIR=$9
TARBALLNAME=${10}
EXECUTALBE_THREADS=${11}
# setup with analysis clone if needed
set -e
source $ANALYSIS_PATH/CROWN/init.sh $ANALYSIS

# use a fourth of the machine for compiling
THREADS_AVAILABLE=$(grep -c ^processor /proc/cpuinfo)
THREADS=$(( THREADS_AVAILABLE / 4 ))
echo "Using $THREADS threads for the compilation"
which cmake

cmake $CROWNFOLDER \
	 -DANALYSIS=$ANALYSIS \
	 -DCONFIG=$CONFIG \
	 -DSAMPLES=$SAMPLES \
	 -DERAS=$ERAS \
	 -DSCOPES=$SCOPE \
	 -DSHIFTS=$SHIFTS \
	 -DTHREADS=$EXECUTALBE_THREADS \
	 -DINSTALLDIR=$INSTALLDIR \
	 -B$BUILDDIR 2>&1 |tee $BUILDDIR/cmake.log

cd $BUILDDIR
make install -j $THREADS
cd $INSTALLDIR
touch $TARBALLNAME
tar -czvf $TARBALLNAME --exclude=$TARBALLNAME .
