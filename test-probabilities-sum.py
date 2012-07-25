import bayes
import math

print("PoissonApprox2 sum over lmbda. Should be ~1 for large k.")
k=1
print("PoissonApprox2 (k=%s): %s" % (k,math.fsum(\
	[bayes.PoissonApprox2(k,lmbda) for lmbda in range(k*3)])))
k=2
print("PoissonApprox2 (k=%s): %s" % (k,math.fsum(\
	[bayes.PoissonApprox2(k,lmbda) for lmbda in range(k*3)])))
k=5
print("PoissonApprox2 (k=%s): %s" % (k,math.fsum(\
	[bayes.PoissonApprox2(k,lmbda) for lmbda in range(k*3)])))
k=10
print("PoissonApprox2 (k=%s): %s" % (k,math.fsum(\
	[bayes.PoissonApprox2(k,lmbda) for lmbda in range(k*3)])))
k=100
print("PoissonApprox2 (k=%s): %s" % (k,math.fsum(\
	[bayes.PoissonApprox2(k,lmbda) for lmbda in range(k*3)])))
k=1000
print("PoissonApprox2 (k=%s): %s" % (k,math.fsum(\
	[bayes.PoissonApprox2(k,lmbda) for lmbda in range(k*3)])))

print("All probabilities that follow should be ~1.")

#P_Nbpi
b=bayes.Bayes()
o=bayes.Observation()
print("P_Nbpi: %s" % math.fsum([b.P_Nbpi(o,Nbpi) for Nbpi in range(3*4000)]))

#P_Npi_S
#to sum this we must produce every sensible number of on counts,
#given that the off counts are 4000 (the default)
#so we'll say max on counts are 4000
#we'll also try two different values of S
b=bayes.Bayes()
b.observations=[bayes.Observation() for i in range(500)]
for i,o in enumerate(b.observations):
	o.Nmi=100
	o.Npi=i
S=b.Slist[3]
print("P_Npi_S: %s" % math.fsum([b.P_Npi_S(o,S) for o in b.observations]))

#P_S_Np_all
#easy to test... just sum the list
b=bayes.Bayes()
o=bayes.Observation()
b.observations=[o]
print("P_S_Np_all: %s" % math.fsum(b.P_S_Np_all()))

