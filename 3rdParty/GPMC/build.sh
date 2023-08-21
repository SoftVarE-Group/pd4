#!/bin/sh

mkdir -p build
mkdir -p bin
cd build

if [ $1 = "d" ]; then
  echo "Create Makefile for Debug..."
  cmake -DCMAKE_BUILD_TYPE=Debug .. \
  && make -j2\
  && cp -p gpmc ../bin/ 

elif [ $1 = "r" ]; then
  echo "Create Makefile for Release..."
  cmake -DCMAKE_BUILD_TYPE=Release .. \
  && make -j2\
  && cp -p gpmc ../bin/ 

elif [ $1 = "rs" ]; then
  echo "Create Makefile for Release Static..."
  cmake -DCMAKE_BUILD_TYPE=Release -DGPMC_STATIC_BUILD=on .. \
  && make -j2\
  && cp -p gpmc ../bin/ 

elif [ $1 = "clean" ]; then
  make clean

fi