from __future__ import print_function # Python2 compatibility
from mcerp import *

x1 = N(24, 1)
x2 = N(37, 4)
x3 = Exp(2)
print(x1.stats)

Z = (x1*x2**2)/(15*(1.5 + x3))
print(Z)
Z.describe()

x1.plot(); x1.show()
Z.plot(); Z.show()
Z.plot(hist=True); Z.show()
Z.plot(); Z.plot(hist=True); Z.show()

print(correlation_matrix([x1, x2, x3]))
plotcorr([x1, x2, x3], labels=['x1', 'x2', 'x3'])
plt.show()

c = np.array([[  1.0, -0.75, 0.0],
              [-0.75,   1.0, 0.0],
              [  0.0,   0.0, 1.0]])

correlate([x1, x2, x3], c)
print(correlation_matrix([x1, x2, x3]))
plotcorr([x1, x2, x3], labels=['x1', 'x2', 'x3'])
plt.show()

Z = (x1*x2**2)/(15*(1.5 + x3))
print(Z)
Z.describe()

print(x1<15)
print(Z>=1000)
