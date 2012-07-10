#python 2.7
import math
from numpy import r_
import scipy.interpolate

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
		kfrac=0.0001
		x=[]
		y=[]
		k=lmbda
		x.append(k)
		p=PoissonApprox2(k,lmbda)
		y.append(p)
		knext=k+1/p*kfrac

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

		print(x)
		print(y)

		self.x=x
		self.y=y

		self.f=scipy.interpolate.interp1d(x,y)
	
	def get(self,k):
		if k>=self.k_max:
			return 0
		if k<0:
			return 0
		return self.f(k)
	
