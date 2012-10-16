import bayes
import time
import math

obs=bayes.Observation()
x=range(4000)
y=[bayes.P_Nnpi(obs, Nnpi) for Nnpi in x]
from pylab import plot, show
plot(x, y)
show()
