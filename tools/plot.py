import matplotlib.pyplot as plt
from numpy import loadtxt, zeros, ones

############################################################################################################################
# PLOTTING PARAMETERS

plt.rcParams[ 'lines.linewidth' ] = 3.5
plt.rcParams[ 'lines.markersize' ] = 12
plt.rcParams[ 'font.size' ] = 28
plt.rcParams[ 'legend.fontsize' ] = 'small'
plt.rcParams[ 'figure.figsize' ] = [ 17, 11 ]
plt.rcParams[ 'axes.labelpad' ] = 18
plt.rcParams[ 'axes.labelsize' ] = 'large'
plt.rcParams[ 'xtick.minor.pad' ] = 15.
plt.rcParams[ 'xtick.major.pad' ] = 15.5
plt.rcParams[ 'ytick.minor.pad' ] = 15.
plt.rcParams[ 'ytick.major.pad' ] = 15.5
plt.rcParams[ 'figure.subplot.bottom' ] = 0.15
plt.rcParams[ 'figure.subplot.left' ] = 0.12
plt.rcParams[ 'figure.subplot.top' ] = 0.96
plt.rcParams[ 'figure.subplot.right' ] = 0.96


############################################################################################################################
# GATHER DATA

zValues    = loadtxt( 'output/zValues.dat' )
analytical = loadtxt( 'output/analytical.dat' )
minVals    = loadtxt( 'output/minVals.dat' )
maxVals    = loadtxt( 'output/maxVals.dat' )




############################################################################################################################
# PLOT

plt.plot( analytical,    zValues, 'k-', label='Analytical' )
plt.plot( minVals,       zValues, 'bo',  label='Min values' )
plt.plot( maxVals,       zValues, 'ro',  label='Max values' )
plt.legend()
plt.grid(True)
plt.xlabel(r'$\mathsf{u}$', fontsize=(plt.rcParams[ 'font.size' ])*1.6)
plt.ylabel(r'$\mathsf{z}$', fontsize=(plt.rcParams[ 'font.size' ])*1.6)
#plt.yscale('log')
#plt.xscale('log')
plt.show()








