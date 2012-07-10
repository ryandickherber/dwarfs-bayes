import bayes
import math

pi=bayes.PoissonInterp(100000)
newx=[k for k in range(300000)]
newy=[pi.get(k) for k in range(300000)]
print("sum: %s" % math.fsum(newy))
print("n points: %s" % len(pi.x))
from pylab import plot, show
plot(newx, newy, pi.x, pi.y, 'o')
show()
