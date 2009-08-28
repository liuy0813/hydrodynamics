##############################################################
#                        Euler1D.py                          #
#                                                            #
#                       19 Feb 2009                          #
#                    by J.C. Toledo-Roy                      #
##############################################################
# Python port of part of the Euler 1D code, implementing the #
# Lax-Friedrichs solver. It was done to compare performance  #
# with the Fortran 90 version.                               #
#                                                            #
# The Sod Shock Tube was run with both codes, with NT=0.5,   #
# CP=0.5, and NX=5000. The final output of both codes was    #
# checked to ensure they produced the same result. The       #
# execution times were:                                      #
#                                                            #
# Python code: 734.75s                                       #
# Fortran90 code: 6.468s                                     #
#                                                            #
# As seen, there is a monumental 100X difference in          #
# performance making Python unsuitable for serious HD work.  #
##############################################################

from __future__ import division
from math import sqrt
from time import strftime, clock
from string import zfill

#######################
# Global Declarations #
#######################

# Define Simulation Parameters
NX = 5000       # Number of grid points
L = 1.0         # Grid physical length
NT = 0.5        # Final integration time
CP = 0.5        # Courant Parameter
GAM = 1.4       # Heat capacity ratios
NOUT = 20       # Number of outputs
BENCH = 1       # Benchmark mode?

# Output Filename Template
# Iteration number will be placed between root and extension
# DIR: data dump directory path (must be created beforehand)
# ROOT: root file name
# EXT: extension
DIR = './data/'
ROOT = 'Sod.'
EXT = '.dat'

# Declare Simulation Variables
# These are lists of lists, structured as follows:
#  U[X][0] = x       P[X][0] = x     F[X][0] = x
#  U[X][1] = rho     P[X][1] = u     F[X][1] = rho*u
#  U[X][2] = rho*u   P[X][2] = rho   F[X][2] = rho*u^2+P
#  U[X][3] = E       P[X][3] = P     F[X][3] = u(E+P)
#  where E = 0.5*rho*u^2 + P/(gamma-1)
U = []     # Integration variables
UP = []    # Stepped integration variables
P = []     # Primitives
F = []     # Fluxes
# ... and initialize them
for X in range(NX):
  U.append([X,0,0,0])
  UP.append([X,0,0,0])
  P.append([X,0,0,0])
  F.append([X,0,0,0])

# Other simulation parameters
IT = 0         # Iteration number
T = 0          # Actual Time
DT = 0         # Timestep size
DX = L / NX    # Grid spacing
NDUMP = 0      # Output number 


#####################
# Program Functions #
#####################

###########################
# Load Initial Conditions #
###########################
def initcond():
   # Sod Shock Tube initial conditions
   # This is a standard test for hydro codes
   # Velocity = 0 everywhere
   # Density = 1.0 on left half, 0.125 on right half
   # Pressure = 1.0 on left, 0.1 on right
  for X in range(NX):
    if (X+1<=NX/2):
      P[X][1] = 0.0
      P[X][2] = 1.0
      P[X][3] = 1.0
    else:
      P[X][1] = 0.0
      P[X][2] = 0.125
      P[X][3] = 0.1
  getvars()

#############################
# Calculate integration     #
# variables from primitives #
#############################
def getvars():
  for X in range(NX):
    U[X][1] = P[X][2]
    U[X][2] = P[X][2]*P[X][1]
    U[X][3] = 0.5*P[X][2]*P[X][1]**2 + P[X][3]/(GAM-1)

##############################
# Calculate primitives from  #
# integration variables      #
##############################
def getprims():
  for X in range(NX):
    P[X][1] = U[X][2]/U[X][1]
    P[X][2] = U[X][1]
    P[X][3] = (GAM-1)*(U[X][3]-0.5*U[X][2]**2/U[X][1]) 

######################
# Calculate Timestep #
######################
def getstep():
  # The timestep is calculated following the
  # Courant-Friedrichs-Lewis criterion:
  # DT = DX / UMAX * CP
  # where DX is the grid spacing, CP is the 
  # Courant parameter and UMAX is the maximum
  # value of local velocity + local sound speed
  UMAX = 0
  ULOC = 0
  for x in range(NX):
    ULOC = abs(P[x][1]) + sqrt(GAM*P[x][3]/P[x][2])
    if (ULOC > UMAX):
      UMAX = ULOC
  return DX/UMAX * CP

##############
# Lax Method #
##############
def lax():

  # Calculate fluxes
  for X in range(NX):
    F[X][1] = P[X][2]*P[X][1]
    F[X][2] = P[X][2]*P[X][1]**2 + P[X][3]
    F[X][3] = P[X][1]*(0.5*P[X][2]*P[X][1]**2 + GAM/(GAM-1)*P[X][3])

  # Calculate timestepped variables
  for X in range(1,NX-1):
    for V in [1,2,3]:
      UP[X][V] =  0.5*(U[X+1][V] + U[X-1][V]) - 0.5 * DT/DX * (F[X+1][V] - F[X-1][V])

  # Treats boundary grid points
  boundary()


#######################
# Boundary Conditions #
#######################
def boundary():

  # Reflection boundary condition
  UP[0][2] = -1*UP[NX-2][2]
  for V in [1,3]:
    UP[0][V] = UP[1][V]
    UP[NX-1][V] = U[NX-2][V]

######################
# Advance Simulation #
######################
def advance():
  global T, IT

  for X in range(NX):
    for V in [1,2,3]:
      U[X][V] = UP[X][V] 
  T = T + DT
  IT = IT + 1

######################
# Dump primitives to #
# output data file   #
######################
def output():
  global NDUMP

  # Build output filename
  fname = DIR+ROOT+zfill(str(NDUMP),4)+EXT
  file = open(fname,"w")

  print strftime('%b %d %Y @ %H:%M:%S') \
    + ' - Dumping Ouput ' + str(NDUMP) \
    + ' to ' + fname

  # Header
  file.write('# Dumped on XXX\n')
  file.write('# Using method XXX\n')
  file.write('# <rest of header>\n')

  # Data
  for X in range(NX):
    s = str(X*DX) + '\t' \
      + str(P[X][1]) + '\t' \
      + str(P[X][2]) + '\t' \
      + str(P[X][3]) + '\n'
    file.write(s)

  # Close file
  file.close()

##########################
# Reports progress and   #
# dumps output if needed #
##########################
def progress():
  global NDUMP

  if (T/NT > NDUMP/NOUT):
    print 'Progress: %6.3f%% (IT=%5.5i)' %  (T/NT*100, IT)
    if (BENCH!=1): output()
    NDUMP = NDUMP + 1

# Benchmarking function
def bench(secs):
  print 'Time elapsed: ' + str(secs) + 's'

#################
# PROGRAM START #
#################

# Initialize
initcond()
if (BENCH!=1): output()

# Timer
if (BENCH==1): print 'Starting simulation in benchmark mode ...'
start = clock()

# Main Loop
while(T<=NT):

  DT = getstep()  # Get timestep

  lax()         # Lax Solver
  advance()     # Advance simulaiton
  getprims()    # Update primitives

  progress()    # Show progress; dump output if needed

end = clock()
bench(end-start)
