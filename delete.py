import numpy as np

i = open('lognormpdf_theory_500p10.txt', 'r')
p = i.read()
pp = p.split()
ppp = np.asarray(pp)
for s in ppp:
    v = complex(s)

print (ppp.shape)
