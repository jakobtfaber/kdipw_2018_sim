
# import statmements for book plotting
# check that these are the best approach
# i.e. should I always use np and eliminate the 'from numpy import *' statement?
# is 'from pylab import *' necessary?
import numpy as np
from math import *
from pylab import *
from scipy import *
from scipy import stats

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import use
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show 

#plt.style.use('seaborn-notebook')	# nice pastel grid
#plt.style.use('classic')		# classic framing
plt.style.use('default')		# default (same as classic?)

# can create custom style sheet: see https://matplotlib.org/users/customizing.html

# use rcParams = dictionary like variable  to control line widths etc.  
# rcParams without args to see long list

# rcdefaults() restores defaults

# rc = convenience function for modifying rc settings. 

# can set global defaults in /Users/cordes/.matplotlib/matplotlibrc
# or can have one in the working directory

# first set defaults

mpl.rcdefaults()
#mpl.rcParams['text.usetex'] = True  	# this seems to be slow

"""
for arrows see:
https://matplotlib.org/users/annotations_guide.html

also:
/opt/local/share/py27-matplotlib/examples/pylab_examples/annotation_demo2.py

note if arrowprops=dict(arrowstyle='<-') is used in the annotate command, then
other arrow properties like headwidth, headlength, width cannot be used. 
"""
