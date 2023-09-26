#!/bin/bash

set -e
set -u
set -o pipefail

opt=0

while getopts 'dsl' OPTION
do
    case "$OPTION" in
        d)
            opt=1
            ;;
        s)
            opt=2
            ;;
        l)
            opt=3
            ;;

    esac
done

curRep=$PWD
echo building glucose
if [ $opt -eq 1 ]
then
   if ! [ -f 3rdParty/glucose-3.0/core/lib_debug.a ]
   then
       cd 3rdParty/glucose-3.0/core/
       make libd       
   fi
elif [ $opt -eq 2 ]
then
    if ! [ -f 3rdParty/glucose-3.0/core/lib_static.a ]
    then
        cd 3rdParty/glucose-3.0/core/
        make libst       
    fi
else
    if ! [ -f 3rdParty/glucose-3.0/core/lib_standard.a ]
    then
        cd 3rdParty/glucose-3.0/core/
        make libs
    fi
fi



echo building kahypar
if ! [ -f 3rdParty/kahypar/build/lib/libkahypar.a ]
then
    cd $curRep
    cd 3rdParty/kahypar/
    mkdir -p build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=RELEASE
    make -j
fi



echo building cryptominisat
if ! [ -f ./3rdParty/cryptominisat/build/lib/libcrytominisat5.so]
then
    cd $curRep
    cd 3rdParty/cryptominisat
    mkdir -p build
    cd build
    cmake -DCMAKE_CXX_FLAGS="-include cstdint"  .. 
    make -j12
fi

echo building louvain-community
if ! [ -f ./3rdParty/louvain-community/build/lib/liblouvain_communities.so  ]
then
    cd $curRep
    cd 3rdParty/louvain-community
    mkdir -p build
    cd build
    cmake -DCMAKE_CXX_FLAGS="-include cstdint"  .. 
    make -j12
fi

echo building arjun
if ! [ -f ./3rdParty/arjun/build/lib/libarjun.so  ]
then
    cd $curRep
    cd 3rdParty/arjun
    mkdir -p build
    cd build
    cmake -DCMAKE_PREFIX_PATH="../../louvain-community/build;../../cryptominisat/build"  ..
    make -j12
fi


echo building mt-kahypar
if ! [ -f 3rdParty/mt-kahypar/build/lib/libmtkahypar.so ]
then
    cd $curRep
    cd 3rdParty/mt-kahypar
    mkdir -p build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=RELEASE  -DKAHYPAR_DOWNLOAD_TBB=On -DKAHYPAR_DOWNLOAD_BOOST=On
    cmake --build . --target mtkahypar
fi

cd $curRep
mkdir -p build
cd build
cmake -GNinja .. -DBUILD_MODE=$opt 
ninja 
