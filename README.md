## StocDeltaN
StocDeltaN is a numerical solver of the stochastic-delta N formalism writen in C++.

## How to use
All required programs are packed into `source` directory. Please download this directory to your working space.

Let's first play with sample codes! Download `sample` directory to the same place of `source` directory.
StocDeltaN is written in C++11, so please use `-std=c++11` option to compile the sample program.  
e.g. for GNU compiler

    g++ -std=c++11 sample_DC_conf.cpp  
    
for Intel compiler

    icpc -std=c++11 sample_DC_conf.cpp  

If you want to hasten a calculation, you can use OpenMP option to parallelize it.  
e.g. for GNU compiler (Note that recent Mac adopts Clang compiler as g++ which does not support OpenMP. To use OpenMP, please manually install GNU compiler or clang-omp.)

    g++ -std=c++11 -fopenmp sample_DC_conf.cpp

for Intel compiler

    icpc -std=c++11 -qopenmp sample_DC_conf.cpp

