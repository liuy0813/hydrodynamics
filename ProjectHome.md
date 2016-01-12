This project contains the programs developed during a seminar that I took in 2009 at ICN-UNAM, titled "Topics in Gasdynamics in Astrophysics". The purpose of the seminar was to expand the basic 1D and 2D HD codes we developed in the first part of the seminar by implementing more advanced solvers as well as parallelization.

# Euler1D Code #

This is a 1D hydrodynamic code which tests several solution methods for the Euler equations, oriented towards Riemann Problems. Visit the [main wiki page](Euler1DWiki.md) for more details.

### Benchmarks ###

In the spirit of comparing programming language performance, I've written ports of the original Fortran90 Euler1D code in C++, Java and Python.

The results were somewhat as expected: Fortran90 was the fastest, with optimized C++ trailing only 3% behind and Java 9% behind (pretty good for a cross-platform, interpreted language!). The Python version, as well as a run of the Java code with the Just-In-Time compiler disabled, were catastrophically slower: 87x times slower than F90 for Java without JIT and 113x for Python. It's quite interesting to have hard numbers on just how much faster Java code is due to the JIT compiler. You can view the details in the [Benchmarks page](Euler1D_Benchmarks.md).

# Euler2D Code #

Development on this "toy code" has stopped to give priority to my actual work projects.