import numpy as np
import matplotlib.pyplot as plt
import math

rand.seed(1234234345)

nrays = 10

x = np.zeros(nrays)
y = np.zeros(nrays)
x2 = np.zeros(nrays)
y2 = np.zeros(nrays)

for i in range(nrays):
    
    r = 1

    theta = rand.gauss(-np.pi/2000, np.pi/2000)
    phi = rand.gauss(-np.pi/2, np.pi/2)
    
    #if -np.pi/60 <= theta <= np.pi/60:

    x[i] = r * np.cos(theta) * np.sin(phi)
    y[i] = r * np.sin(theta) * np.sin(phi)
    
    for i in range(nrays):
        
        theta2 = rand.gauss(-np.pi/2, np.pi/2)
        phi2 = rand.gauss(-np.pi/2, np.pi/2)
        
        x2[i] = ((r) * np.cos(theta2) * np.sin(phi2)) + x[i]
        y2[i] = ((r) * np.sin(theta2) * np.sin(phi2)) + y[i]

            
            
fig = plt.figure(figsize = (12,8))
ax1 = fig.add_subplot(121)
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.scatter(y,x, c = 'black', marker = 'o')
plt.gca().set_aspect('equal')

ax2 = fig.add_subplot(122)
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.scatter(y2,x2, c = 'black', marker = 'o')
plt.gca().set_aspect('equal')



plt.show()