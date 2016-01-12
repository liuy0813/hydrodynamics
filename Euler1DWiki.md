### Contents ###



# Introduction #

Euler1D is a rudimentary 1-dimensional hydrodynamics code, originally written as part of a Seminar on Numerical Gas Dynamics that I took in August-December 2007, and expanded during a follow-up seminar in early 2009.

The code solves the 1D inviscid Euler equations in conservation form by means of a variety of first- and second-order methods. It was tested with a small battery of shock-inducing problems given by E.F. Toro ([book](http://www.amazon.com/Riemann-Solvers-Numerical-Methods-Dynamics/dp/3540616764/ref=sr_1_6?ie=UTF8&s=books&qid=1238732924&sr=1-6)), which are similar to _Sod's Shock Tube_.

The code is far from being complete, lacking many features that would make it more portable (such as better handling of numerical divergence, creating the data/ dir automatically, etc). However, it should still be able to simulate any decently conceived 1D problem.

# Details #

## The Physics ##

The 1D Euler equations for inviscid flows in vector conservation form are:

<img src='http://img6.imageshack.us/img6/5681/latex2png2.png'>

where<br>
<br>
<img src='http://img6.imageshack.us/img6/3244/latex2png2h.png'>

Here, <img src='http://img154.imageshack.us/img154/5681/latex2png2.png'> is the flow velocity, <img src='http://img154.imageshack.us/img154/9494/latex2png2w.png'> is the fluid mass density, <img src='http://img518.imageshack.us/img518/5681/latex2png2.png'> is the gas pressure and<br>
<br>
<img src='http://img11.imageshack.us/img11/5681/latex2png2.png'>

is the total energy per unit volume (with <img src='http://img48.imageshack.us/img48/5681/latex2png2.png'> the internal thermal energy per unit volume). The system of equations is not closed and requires an equation of state to be solvable. In this case we study only ideal gasses, with an equation of state <img src='http://img237.imageshack.us/img237/5681/latex2png2.png'>. The total energy then reduces to<br>
<br>
<img src='http://img141.imageshack.us/img141/5681/latex2png2.png'>

where <img src='http://img21.imageshack.us/img21/5681/latex2png2.png'> is the heat capacity ratio. All simulations on this page use 1.4 for this parameter.<br>
<br>
<h2>Solution Algorithms</h2>

The code implements the following algorithms to solve the Euler equations:<br>
<br>
<ul><li>Finite-difference methods<br>
<ol><li>Lax-Friedrichs method (<a href='http://en.wikipedia.org/wiki/Lax%E2%80%93Friedrichs_method'>Wikipedia</a>)<br>
</li><li>MacCormack method (<a href='http://en.wikipedia.org/wiki/MacCormack_method'>Wikipedia</a>)<br>
</li><li>MacCormack method with flux-corrected viscosity<br>
</li></ol></li><li>Finite-element methods<br>
<ol><li>Simple Godunov upwind scheme<br>
</li><li>HLL Riemann solver (Harten, Lax, van Leer (1983))<br>
</li><li>HLLC Riemann solver<br>
</li><li>HLLC with Linear Data Reconstruction (experimental)</li></ol></li></ul>

In addition, the HLL/HLLC solvers can use a variety of wavespeed estimation methods:<br>
<br>
<ul><li>Four simple <i>direct</i> estimators<br>
</li><li>Three <i>pressure-based</i> estimators<br>
</li><li>An adaptive non-iterative pressure-based estimator</li></ul>

<h2>Tests</h2>

The code has been subject to 5 tests proposed by E.F. Toro. In all tests, the physical domain (of length 1.0) is initialized with different, but uniform, values on each side of a central discontinuity. The specific values for density, pressure and velocity are given in the following table.<br>
<br>
<table><thead><th> <b>Test</b> </th><th> <img src='http://img516.imageshack.us/img516/5681/latex2png2.png'> </th><th> <img src='http://img516.imageshack.us/img516/3649/latex2png2r.png'> </th><th> <img src='http://img520.imageshack.us/img520/5681/latex2png2.png'> </th><th> <img src='http://img520.imageshack.us/img520/5722/latex2png2e.png'> </th><th> <img src='http://img520.imageshack.us/img520/1525/latex2png2a.png'> </th><th> <img src='http://img128.imageshack.us/img128/5681/latex2png2.png'> </th></thead><tbody>
<tr><td> 1           </td><td> 1.0                                                                </td><td> 0.75                                                                </td><td> 1.0                                                                </td><td> 0.125                                                               </td><td> 0.0                                                                 </td><td> 0.1                                                                </td></tr>
<tr><td> 2           </td><td> 1.0                                                                </td><td> -2.0                                                                </td><td> 0.4                                                                </td><td> 1.0                                                                 </td><td> 2.0                                                                 </td><td> 0.4                                                                </td></tr>
<tr><td> 3           </td><td> 1.0                                                                </td><td> 0.0                                                                 </td><td> 1000.0                                                             </td><td> 1.0                                                                 </td><td> 0.0                                                                 </td><td> 0.01                                                               </td></tr>
<tr><td> 4           </td><td> 5.99924                                                            </td><td> 19.5975                                                             </td><td> 460.894                                                            </td><td> 5.99242                                                             </td><td> -6.19633                                                            </td><td> 46.0950                                                            </td></tr>
<tr><td> 5           </td><td> 1.0                                                                </td><td> -19.59745                                                           </td><td> 1000.0                                                             </td><td> 1.0                                                                 </td><td> -19.59745                                                           </td><td> 0.01                                                               </td></tr></tbody></table>

<br>

<h1>Installation and Usage</h1>

<h2>Requirements</h2>

To compile and run the code you will need:<br>
<ul><li>A Fortran 90 compiler. I used the GCC Fortran 4.2.1 compiler (gfortran).<br>
</li><li>The PGPLOT Graphics Subroutine Library, available <a href='http://www.astro.caltech.edu/~tjp/pgplot/'>here</a>. The tests were run with PGPLOT 5.2 installed.<br>
</li><li>Alternatively, you may wish to comment out all the lines containg PGPLOT subroutine calls (they begin with PG) and the call to the PLOT subroutine (sorry for the lack of preprocessor directives to handle this), and plot the data manually. The code dumps the simulation data at regular intervals to the data/ directory, in a format suitable to be read and plotted by gnuplot (and I guess other programs).</li></ul>

<h2>Installing</h2>

Installing and running the Euler1D code is simple:<br>
<ul><li>Download the source code from <a href='http://code.google.com/p/hydrodynamics/source/browse/trunk/Euler1D/'>here</a>.<br>
</li><li>Also, download the tictoc.f90 source from <a href='http://code.google.com/p/hydrodynamics/source/browse/#svn/trunk/lib'>here</a>. Put this file in the same directory (it's a small library to time the execution of the code).<br>
</li><li>Compile the code by running <code>make</code>. If everything went well, you should have an 'Euler1D' binary.<br>
</li><li><font color='red'>Important:</font> create a directory named 'data' before running the code, or you will get an error.<br>
</li><li>Run the code by executing the binary: <code>./Euler1D</code></li></ul>

You should get a lot of feedback in the terminal, indicating progress and data dumps, as well as a real-time plot in a separate PGPLOT window (you will also notice a small black window, possibly in the upper left corner: this is the PGPLOT server). The simulation data will be written to separate files in the data/ directory.<br>
<br>
The test run executed by default (as of <a href='https://code.google.com/p/hydrodynamics/source/detail?r=9'>r9</a>) is Test 1 of Toro, which is very similar to Sod's Shock Tube: an initial discontinuity in both pressure and density which evolves into a right-moving shock wave and a left-moving rarefaction wave. Here's a sample of how it should look in the PGPLOT window:<br>
<br>
<a href='http://img123.imageshack.us/img123/2838/screenyvt.jpg'>

<IMG src=http://img123.imageshack.us/img123/2838/screenyvt.th.jpg border="0">

<a />

<h2>Settings</h2>

If you wish to change the settings of the code, or choose a different test, you must modify the <b>params.f90</b> source file, and recompile the code. This configuration file is heavily commented and should be self-explanatory. Amongst the things you can change are:<br>
<br>
<ul><li>The number of grid points and size of the grid.<br>
</li><li>The Courant parameter of the simulation.<br>
</li><li>The solver used to integrate the Euler equations.<br>
</li><li>The wavespeed estimation method for HLL/HLLC solvers.<br>
</li><li>The type of boundary conditions: transmission, reflection, periodic.<br>
</li><li>The initial conditions. The code has 5 built-in tests (ICS=1 to 5).<br>
</li><li>The destination and filename template for data dumps.