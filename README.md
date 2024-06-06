# pd4 - An Extension of d4 for Projected d-DNNF Compilation

# How to Compile


## Dependencies
This repository includes several sub-repositories. So, make sure to clone the repository with:
```console
$ git clone --recurse-submodules <GIT-URL>
```
In addition, all of the listed dependencies are required to compile pd4.

* Build system
  * cmake (Version >= 3.1)
  * ninja
  * gcc/g++
* Libraries
  * zlib
  * m4
  * gmp
  * boost (static)
  * hwloc
  * lzma
  * bzip2
  * libstd
  * tbb

## Building
To build the project execute the build script which also handles building the submodules.
```console
$ ./build.sh
```

The executable is called d4 and is in the build repository.

```console
$ ./build/d4 -h
```

## Running

### Input Format
The input format for projected d-DNNF compilation is DIMACS with an established extension for projected variables.
With the extension, the second line indicates the variable to keep during the projection. Invoking pd4 with the listed input will provide a projected d-DNNF only with the variables 1 2 3.

```
p cnf 4 3
c p show 1 2 3 0
3 1 2 0
−2 1 3 0
2 −3 2 4 0
```

### Command
This command produces a projected d-DNNF for the input file `input.cnf`
```console
./build/d4 -i input.cnf -m proj-ddnnf-compiler 
```