#python 2.7
import pickle
import time
import math
from numpy import r_
import scipy.interpolate

class Observation:
	def __init__(self):
		#TODO: fill in realistic default values
		self.Ci=4.0
		self.Ti=20*60 #seconds
		self.Npi=1000 #on counts
		self.Nmi=4000 #off counts
		self.Ri=150.0 #reconstructed energy average
		#self.Ri0=100.0 #reconstructed energy bottom GeV
		#self.Ri1=200.0 #reconstructed energy top GeV
		self.Ei=[50.0, 150.0, 250.0] #true energy average in GeV
		self.Eirange=[(10.0,100.0),(100.0,200.0),(200.0,400.0)]
		self.P_Ei=[0.25, 0.5, 0.25]
		#avg effective areas in each range, meters squared
		self.Ai=[\
			[5000.0],\
			[10000.0],\
			[11000.0],\
			]
		self.Airange=[\
			[(4000.0,6000.0)],\
			[(9000.0,11000.0)],\
			[(10000.0,12000.0)],\
			]
		self.P_Ai=[\
			[1.0],\
			[1.0],\
			[1.0],\
			]
		#J factors, measured in whatever units they are usually measured in
		self.Ji=[5, 6, 7]
		self.Jirange=[(4.5,5.5),(5.5,6.5),(6.5,7.5)]
		self.P_Ji=[0.25, 0.5, 0.25]

class Bayes:
	def __init__(self, Mchi=400.0, Slist=[0,1,2,3,4,5,6,7,8,9,10]):
		self.Mchi=Mchi #GeV #TODO: fill in realistic value for mass

		#TODO: fill in realistic values for cross section
		self.Slist=Slist

		self.observations=[]
		self.P_S_sum=0
		self.memoizedict={}

	def P_Nnpi(self, obs, Nnpi):
		prob=PoissonApprox2(obs.Nmi,obs.Ci*Nnpi)
		#print("P_Nnpi: %s" % prob)
		return prob

	def P_Npi_JiAiEiNnpiS(self, obs, iJi, iAi, iEi, Nnpi, S):
		Mchi=self.Mchi
		Ngi=S/(8*math.pi*(Mchi**2))*\
			integrate_dNdE(obs.Eirange[iEi][0],obs.Eirange[iEi][1])\
			*obs.Ai[iEi][iAi]*obs.Ti*obs.Ji[iJi]
		lmbda=Nnpi+Ngi
		k=obs.Npi
		prob=PoissonApprox2(k,lmbda)
		#print("P_Npi_JiAiEiNnpiS: %s" % prob)
		return PoissonApprox2(k,lmbda)

	def P_Npi_NnpiS(self, obs, Nnpi, S):
		p=[]
		for iEi, Ei in enumerate(obs.Ei):
			for iAi, Ai in enumerate(obs.Ai[iEi]):
				for iJi, Ji in enumerate(obs.Ji):
					P=self.P_Npi_JiAiEiNnpiS(obs, iJi, iAi, iEi, Nnpi, S)
					P=P*obs.P_Ji[iJi]*obs.P_Ai[iEi][iAi]*obs.P_Ei[iEi]
					p.append(P)
		prob=math.fsum(p)
		#print("P_Npi_NnpiS: %s" % prob)
		return math.fsum(p)

	def P_Npi_S(self, obs, S):
		Nnpi_cutoff=int(3*obs.Nmi/obs.Ci)
		if Nnpi_cutoff<30:
			Nnpi_cutoff=30
		p=[self.P_Nnpi(obs, Nnpi)*self.P_Npi_NnpiS(obs, Nnpi, S)\
			for Nnpi in range(Nnpi_cutoff)]
		prob=math.fsum(p)
		#print("P_Npi_S: %s" % prob)
		return prob

	def P_S_Np_product(self, S):
		Mchi=self.Mchi
		T=[]
		for obs in self.observations:
			P=self.P_Npi_S(obs, S)
			if P!=0:
				T.append(math.log(P))
		t=math.fsum(T)
		t=math.exp(t)
		return t

	def P_S_Np_product_memoized(self, S):
		Mchi=self.Mchi
		if S in self.memoizedict:
			if Mchi in self.memoizedict[S]:
				return self.memoizedict[S][Mchi]
			else:
				self.memoizedict[S][Mchi]=self.P_S_Np_product(S)
				return self.memoizedict[S][Mchi]
		self.memoizedict[S]={}
		self.memoizedict[S][Mchi]=self.P_S_Np_product(S)
		return self.memoizedict[S][Mchi]

	def P_S_Np(self, S):
		Mchi=self.Mchi
		top=self.P_S_Np_product_memoized(S)
		bottom=0
		B=[]
		if self.P_S_sum==0:
			for S in self.Slist:
				t=self.P_S_Np_product_memoized(S)
				B.append(t)
			self.P_S_sum=math.fsum(B)
		bottom=self.P_S_sum
		return top/bottom

def integrate_dNdE(Ei0, Ei1):
	#TODO: make this work
	return 1.0

def Poisson(k,lmbda):
	if (lmbda==0):
		return 0
	return (lmbda**k * math.exp(-lmbda))/(math.factorial(k))

def PoissonApprox1(k,lmbda):
	#first order stirling approxmation... not good enough
	if (lmbda==0):
		return 0
	if (k==0):
		return math.exp(-lmbda)
	stirling=k*math.log(k)-k
	return math.exp(k*math.log(lmbda)-lmbda-stirling)

def PoissonApprox2(k,lmbda):
	#higher order stirling approximation... good enough
	#this is the function that should be used for a Poisson
	if (lmbda==0):
		return 0
	if (k==0):
		return math.exp(-lmbda)
	stirling=k*math.log(k)-k+0.5*math.log(2*math.pi*k) \
		+1.0/(12*k)-1.0/(360*(k**3))
	return math.exp(k*math.log(lmbda)-lmbda-stirling)

def PoissonBetterApprox(k,lmbda):
	#higher order stirling approx plus real poisson and low values
	#note that this function is much slower and not more accurate
	#than PoissonApprox2
	if (lmbda**k>10**200 or math.factorial(k)>10**200):
		return PoissonApprox2(k,lmbda)
	else:
		return Poisson(k,lmbda)

class PoissonInterp:
	def __init__(self,lmbda):
		k_max=lmbda*3
		self.k_max=k_max
		kfrac=0.01
		x=[]
		y=[]
		k=lmbda
		x.append(k)
		p=PoissonApprox2(k,lmbda)
		y.append(p)
		knext=k+1/p*kfrac

		t1=time.clock()
		#right side
		while knext<k_max:
			if p!=0:
				knext=int(math.floor(k+1/p*kfrac))
			else:
				knext=k_max
			if k==knext:
				knext=k+1
			if knext>=lmbda*3:
				if x[-1]!=k_max:
					knext=k_max
			x.append(knext)
			pnext=PoissonApprox2(knext,lmbda)
			y.append(pnext)
			k=knext
			p=pnext
		t2=time.clock()
		print("right side time: %ss" % (t2-t1,))

		t1=time.clock()
		#left side
		k=lmbda
		p=PoissonApprox2(k,lmbda)
		while knext>0:
			if p!=0:
				knext=int(math.ceil(k-1/p*kfrac))
			else:
				knext=0
			if k==knext:
				knext=k-1
			if knext<0:
				if x[0]==0:
					break
				else:
					knext=0
			x.insert(0,knext)
			pnext=PoissonApprox2(knext,lmbda)
			y.insert(0,pnext)
			k=knext
			p=pnext
		t2=time.clock()
		print("left side time: %ss" % (t2-t1,))

		#print(x)
		#print(y)

		self.x=x
		self.y=y

		t1=time.clock()
		self.f=scipy.interpolate.interp1d(x,y)
		t2=time.clock()
		print("set up f time: %ss" % (t2-t1,))
	
	def get(self,k):
		if k>=self.k_max:
			return 0
		if k<0:
			return 0
		return self.f(k)

if __name__=="__main__":
	b=Bayes()
	b.observations=[Observation()]
	for S in b.Slist:
		P=b.P_S_Np(S)
		print("S, P: %s, %s" % (S,P))
	print("Sum over P: %s" % math.fsum(\
		[b.P_S_Np(S) for S in b.Slist]))
