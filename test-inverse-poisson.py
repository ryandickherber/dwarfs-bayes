import bayes
import math
import time

k=1000
if (k>10):
    lmbda_max=k*3
else:
    lmbda_max=30
x=range(lmbda_max)
y=[bayes.PoissonApprox2(k,lmbda) for lmbda in x]
from pylab import plot, show
plot(x, y)
show()
