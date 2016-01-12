## Python port of Euler1D ##

Since I love Python, I decided to write a simple Python port of the core Euler1D code, implementing only the Lax-Friedrichs algorithm. The purpose of doing this, of course, was to compare the performance of the port with the original Fortran90 version. I was expecting Python to run considerably slower, but the results were far worse.

You can download the Python port [here](http://code.google.com/p/hydrodynamics/source/browse/trunk/Euler1D/python/Euler1D.py).

### Benchmark ###

I ran the standard [Sod Shock Tube](http://en.wikipedia.org/wiki/Sod_Shock_Tube) test with both codes, using the Lax-Friedrichs algorithm and the following simulation parameters:

  * Number of grid points: 1000
  * Final integration time: 5.0
  * Courant parameter: 0.5

Reflection boundary conditions were used to keep the problem 'alive' up to t=5.0.

Here are the execution times for the test:

  * Python code: 235.72 s
  * Fortran90 code: 2.315 s

I think these numbers speak for themselves. This definitely makes Python unsuitable for serious hydrodynamic simulations.

This is mainly due to the fact that Python is an interpreted (as opposed to compiled) language. The difference is that, when a program is run, the Python interpreter has to parse and check each line of code every time it is executed. This small overhead quickly adds up when you have large loops of instructions, which is typical in hydro codes. With compiled languages, all the code is parsed and checked only once.

Another possible reason Python is slower is because several of the higher-level abstract constructs supported in Python (that make coding a pleasure) are probably slower than the lower-level equivalent implementation in Fortran.

Now, fostering my undiminished love for Python here, I'll say a few things in its defense. First, of course, the comparison I'm making here is very unfair, since Python code is not compiled. It's like comparing a limousine to a sport car - they have different design goals. While there are small 3rd party Python compilers, I don't think they can compete with the GCC compiler: it is heavily optimized, having been improved continuously over decades by a large base of developers.

Second, I'm sure my Python implementation of the Euler1D code could be optimized (I hacked it together in very little time - that's the strong point of Python). For instance, I use the `append` list method a lot, which is probably slower than having static arrays (i.e. lists). However, I suspect that no matter how much optimization is done on the code, it will never approach Fortran's performance (and I'm sure the C performance would be as good).

So Python for fun, Fortran/C for work.