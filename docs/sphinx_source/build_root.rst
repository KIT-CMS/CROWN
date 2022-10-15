How to build ROOT with CVMFS on CentOS 7
=========================================

For profiling with debug symbols or just to test the newest ROOT features, you may want to use your own ROOT version. Here are the commands, which allow you to build ROOT with a given build type, ROOT release tag and C++ 17.

Most likely, you want to use :code:`RelWithDebInfo` as build type so you get debug symbols but also a realistic performance due to compiler optimizations.

To look up the release tags, go to https://github.com/root-project/root and see the tags (not the branches!).


.. code-block:: console

    # Select the release tag (the tag, not the branch!) and clone ROOT
    RELEASE_TAG=v6-24-00
    git clone https://github.com/root-project/root -b $RELEASE_TAG --depth 1

    # Source a proper compiler from CVMFS, which supports C++17
    source /cvmfs/sft.cern.ch/lcg/releases/gcc/10.1.0/x86_64-centos7/setup.sh

    # Build ROOT with the given build type and without tmva and roofit
    BUILD_TYPE=RelWithDebInfo
    mkdir build
    cd build
    cmake3 ../root -Dtmva=OFF -Droofit=OFF -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_CXX_STANDARD=17

    # Build!
    make -j96 # 96 cores? Sure!

    # NOTE: Remember, you have to source the compiler and build/bin/thisroot.sh to get everything up and running!
