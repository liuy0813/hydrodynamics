## Euler1D Language Benchmarks ##

The Euler1D code was a good opportunity to test the performance of different programming languages for scientific computing. A small portion of the original Euler1D code was ported first to Python and later to Java and C++. The same problem was solved on all versions of the code, with the results checked for consistency. The results are presented in this page.

### C++, Java and Python ports ###

Specifically, the ports are a trimmed-down version of the Euler1D code featuring only one solution algorithm: the Lax-Friedrichs method. Since this method is simple, porting the code was relatively quick (more so for Python than for Java or C++).

You can download the Python port [here](http://code.google.com/p/hydrodynamics/source/browse/trunk/benchmarks/Euler1D.py).

The Java port is composed of the [Euler1D.java](http://code.google.com/p/hydrodynamics/source/browse/trunk/benchmarks/Euler1D.java) class; it must be launched as main.Euler1D since it's part of the main package.

The C++ port can be obtained [here](http://code.google.com/p/hydrodynamics/source/browse/trunk/benchmarks/Euler1D.c%2B%2B).


### Benchmark ###

I ran the standard [Sod Shock Tube](http://en.wikipedia.org/wiki/Sod_Shock_Tube) test with all three codes, using the Lax-Friedrichs algorithm and the following simulation parameters:

  * Number of grid points: 5000
  * Final integration time: 0.5
  * Courant parameter: 0.5

### Results ###

The following table summarises the various runs, sorted by increasing execution time.

| **Language**  | **Compiler/Interpreter** | **Options** | **Execution Time** | **Ratio to fastest** |
|:--------------|:-------------------------|:------------|:-------------------|:---------------------|
| Fortran90     | GCC Fortran compiler 4.3.3 | -O3         | 6.634s             | 1.0                  |
| C++           | GCC C++ compiler 4.3.3   | -O3         | 6.841s             | 1.031                |
| Java          | Sun JRE 1.6.0\_14        | -O          | 7.280s             | 1.097                |
| C++           | GCC C++ compiler 4.3.3   | None        | 10.606             | 1.5987               |
| Java - no JIT | Sun JRE 1.6.0\_14        | -Djava.compiler=NONE | 569.43s            | 87.2                 |
| Python        | Python 2.6.2             | -O          | 751.9s             | 113.3                |

Note: for Fortran90, C++ and Java, the execution time shown is the average of 5 consecutive runs (for Java without JIT and Python there was no point).

### Discussion ###