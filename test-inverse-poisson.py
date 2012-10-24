import bayes
import math
import time

k=1000
print("k=1000")
if (k>10):
    lmbda_max=k*3
else:
    lmbda_max=30
x=range(lmbda_max)
y=[bayes.PoissonApprox2(k,lmbda) for lmbda in x]
print("Sum: %s" % (math.fsum(y)))
from pylab import plot, show
plot(x, y)
show()

k=5
print("k=5")
if (k>10):
    lmbda_max=k*3
else:
    lmbda_max=30
x=range(lmbda_max)
y=[bayes.PoissonApprox2(k,lmbda) for lmbda in x]
print("Sum: %s" % (math.fsum(y)))
from pylab import plot, show
plot(x, y)
show()

k=2
print("k=2")
if (k>10):
    lmbda_max=k*3
else:
    lmbda_max=30
x=range(lmbda_max)
y=[bayes.PoissonApprox2(k,lmbda) for lmbda in x]
print("Sum: %s" % (math.fsum(y)))
from pylab import plot, show
plot(x, y)
show()

k=1
print("k=1")
if (k>10):
    lmbda_max=k*3
else:
    lmbda_max=30
x=range(lmbda_max)
y=[bayes.PoissonApprox2(k,lmbda) for lmbda in x]
print("Sum: %s" % (math.fsum(y)))
from pylab import plot, show
plot(x, y)
show()

k=0
print("k=0")
if (k>10):
    lmbda_max=k*3
else:
    lmbda_max=30
x=range(lmbda_max)
y=[bayes.PoissonApprox2(k,lmbda) for lmbda in x]
print("Sum: %s" % (math.fsum(y)))
from pylab import plot, show
plot(x, y)
show()

