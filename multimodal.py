import numpy as np
import matplotlib.pyplot as plt
from pylab import linspace

def f(t):
	return 1 / (1.1 - np.exp(-2*t) * np.sin(5*np.pi*t) * np.sin(5*np.pi*t))

t1 = linspace(0.0, 1.0, 200.0)

min1 = 0.2 
min2 = 0.4
min3 = 0.6
min4 = 0.8

plt.plot(t1, f(t1), 'k')

plt.ylabel('f(x)', fontsize=12)
plt.xlabel('x', fontsize=12)

plt.plot(min1, f(min1),'r',marker='o',clip_on=False)
plt.plot(min2, f(min2),'r',marker='o',clip_on=False)
plt.plot(min3, f(min3),'r',marker='o',clip_on=False)
plt.plot(min4, f(min4),'r',marker='o',clip_on=False)

# plt.title("", fontsize = 12)
plt.title("Multimodal Problem", fontsize=12)

plt.axis([0.05, 0.9, 0.0, 4.0])

plt.show()