Michal's Masters Project
=======

These are the sources for my Masters Project in Physics @ University of Edinburgh. As of the time of writing, the project is a clever integrator for a nonlinear numerical equation, for more info see the paper.

## C++ part

The C++ part uses cmake, so to configure and compile use (from inside of C directory)

    mkdir build ; cd build ; cmake .. ; make -j3

### Tests

Tests use the [Boost Test Framefork](http://www.boost.org/doc/libs/1_54_0/libs/test/doc/html/index.html).

## TODO

[ ] Tests work and test the recent code version (including Searcher)
[ ] Test computer instantiation
[ ] Test each computer in a basic and more complex ways
[ ] Check if command line parameters setting works
[ ] Test if searcher works correctly with basic functions

[ ] Make searcher register its parameters in global params
[ ] A bug causing Searcher to diverge is found
[ ] Add more comments + project page
[x] Searcher works
[x] Multiple threads at instantiation with OpenMP
[ ] Fix the run condition

