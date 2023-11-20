
**Please consider using the Python implementation PyCGTOOL at [https://github.com/jag1g13/pycgtool](https://github.com/jag1g13/pycgtool) as it has several advantages over this version.**

# CGTOOL v0.5
## About the project
The aim of this project is to provide a tool to aid in parametrising coarse-grained (CG) molecular mechanics models.  CGTOOL generates coarse-grained models from atomistic trajectories using a user-provided mapping.  Equilibrium values and force constants of bonded terms are calculated by Boltzmann Inversion of histograms collected from the input trajectory allowing good replication of target properties.
Input is GROMACS XTC and GRO files, along with a custom config file based which specifies the mapping and other parameters.
The output is a GROMACS/MARTINI style forcefield file and a coarse-grained representation of the system created in GROMACS GRO format to allow initial simulations to be started with minimal hand tuning.

Experimental support is present for output in LAMMPS formats to be used with the ELBA forcefield, however this is not yet fully implemented.

These programs are work-in-progress.

Input is GROMACS XTC and GRO files, along with a custom config file based which specifies options for the analysis.  Membrane properties can be averaged over the complete simulation trajectory or output in batches.  An example config file lists the available options.

## Dependencies
To download and build CGTOOL, the following is required
* git installed on the system
* CMake 3.12 or newer
* C++ compiler supporting the C++11 standard
* Boost C++ libraries (program\_options module recommended but optional)
* Optional: Doxygen to build developer documentation

## Install
The following works on all UNIX-like systems (Linux / macOS) 
* Clone this repository
```
git clone https://github.com/qiUip/cgtool.git
```
To compile the program:
* `cd` into the main distribution directory 'cgtool', make a new directory 'build' and `cd` into it
```
cd cgtool && mkdir build && cd build
```
* From the build directory execute `cmake ..` or alternatively use the terminal-based GUI with `ccmake ..`
* To build the executable and library execute `make`.  After compilation the executable can be found in the 'bin' directory and the library in the 'lib' directory.
* To build the developer documentation execute `make doc` (requires Doxygen)

## Use
### CGTOOL
* The program requires a CFG config file, a XTC file and a GROMACS GRO file. To call it use the following command (order of options not important):
```
cgtool -c <CFG file> -x <XTC file> -g <GRO file>
```
* An optional GROMACS ITP file may be provided with the `-i <itp file>` option to allow calculation of charges
* The config file specifies the mapping to be applied, an example is present in the test\_data directory
* Help text is available with `cgtool -h` or `cgtool --help`

### RAMSi
* The program requires a CFG config file, a XTC file and a GROMACS GRO file. To call it use the following command (order of options not important):
```
ramsi  -c <CFG file> -x <XTC file> -g <GRO file>
```
* The config file specifies the mapping to be applied, an example is present in the test\_data directory
* Help text is available with `ramsi -h` or `ramsi --help`

### Testing
Currently only a few source files have complete unit tests.  Integration tests ensure that the program compiles successfully and produces correct output for a test dataset.

An option to test the build upon compilation is available with `make check`. Alternatively, to test after compilation execute `ctest`, or `ctest -V` for more verbose testing.

### Contact
* James Graham: <J.A.Graham@soton.ac.uk>
