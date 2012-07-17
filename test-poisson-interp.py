import bayes
import math
import time

lmbda=100
k_max=300
x=[k for k in range(k_max)]

print("Uninterpolated...")
t1=time.clock()
y=[bayes.PoissonApprox2(k,lmbda) for k in range(k_max)]
ysum=math.fsum(y)
t2=time.clock()
print("time: %ss" % (t2-t1,))
print("sum: %s" % ysum)
print("n poisson: %s" % len(x))
print("Interpolated...")
t1=time.clock()
pi=bayes.PoissonInterp(lmbda)
y=[pi.get(k) for k in range(k_max)]
t2=time.clock()
print("time: %ss" % (t2-t1,))
print("sum: %s" % math.fsum(y))
print("n points: %s" % len(pi.x))
from pylab import plot, show
plot(x, y, pi.x, pi.y, 'o')
show()
