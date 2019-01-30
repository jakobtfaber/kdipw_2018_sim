import numpy as np
import matplotlib.pyplot as plt

#Read in xxx text file (dependent on number of realizations)

g = open('xxx_250p10.txt','r')
x = g.read()
xx = x.split()
xxx = np.asarray(xx, dtype = 'f')

"""#Read in xpdf file (dependent on number of realizations)

f = open('xpdf_500p10.txt','r')
l = f.read()
ll = l.split()
lll = np.asarray(ll, dtype = 'f')

#Make array for numerical negative of lll

negl = np.exp(np.negative(lll))

#Read in lognormpdf_theory file (dependent on number of realizations)

i = open('lognormpdf_theory_500p10.txt','r')
p = i.read()
pp = p.split()
ppp = np.asarray(pp)
# Read array as complex numbers

pppp = map(complex, ppp)
"""
print(xxx)

#Plot histogram, exp line, and lognormpdf_theory lines repsectively

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.title("PhiF_d = 10 ::: 250 Realizations")
ax1.set_xlabel("Intensity")
ax1.set_ylabel("PDF f_I(I)")
#ax1.set_xlim([0, 0.5])
ax1.set_yscale("log", nonposy='clip')
#ax1.set_ylim([10**-10, 10**8])
#plt.yscale('log')
plt.hist(xxx, bins=500)
#plt.plot(lll, negl, c = 'm')
#plt.plot(lll, ppp, c = 'g')

#lin-lin representation of histogram (not necessary)

"""fig = plt.figure()
ax2 = fig.add_subplot(111)
plt.title("PhiF_d = 1.0 ::: 60 Realizations")
ax2.set_xlabel("Intensity")
ax2.set_ylabel("PDF f_I(I)")
ax2.set_ylim(0, 500000)
plt.hist(xxx, bins = 500)"""

plt.show()
