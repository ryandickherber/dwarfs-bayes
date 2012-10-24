#python 2.7
import argparse
import pickle
import time
import math
from numpy import r_
import scipy.interpolate
import constants

class Observation:
    def __init__(self):
        self.RunNum=0
        self.EnergyBinIndex=0
        self.EnergyBin=[200.0,200.0+200.0*0.15]

        #TODO: fill in realistic default values
        self.Ci=4.0
        self.Ti=20*60 #seconds
        self.Npi=1000 #on counts
        self.Nmi=4000 #off counts
        #self.Ri=150.0 #reconstructed energy average
        #self.Ri0=100.0 #reconstructed energy bottom GeV
        #self.Ri1=200.0 #reconstructed energy top GeV
        self.EnergyBinIndexes=[0,1,2]
        self.Ei=[50.0, 150.0, 250.0] #true energy average in GeV
        self.Eirange=[(10.0,100.0),(100.0,200.0),(200.0,400.0)]
        self.P_Ei=[0.25, 0.5, 0.25]
        #avg effective areas in each range, meters squared
        self.Ai=[\
            [5000.0],\
            [10000.0],\
            [11000.0],\
            ]
        #self.Airange=[\
        #    [(4000.0,6000.0)],\
        #    [(9000.0,11000.0)],\
        #    [(10000.0,12000.0)],\
        #    ]
        self.P_Ai=[\
            [1.0],\
            [1.0],\
            [1.0],\
            ]
        #J factors, measured in whatever units they are usually measured in
        self.Ji=[5, 6, 7]
        self.Jirange=[(4.5,5.5),(5.5,6.5),(6.5,7.5)]
        self.P_Ji=[0.25, 0.5, 0.25]

    def save(self, filename):
        f=open(filename, 'wb')
        pickle.dump(self, f)
        f.close()
    
    def open(self, filename):
        f=open(filename, 'rb')
        o=pickle.load(f)
        f.close()
        return o

class Bayes:
    def __init__(self, Mchi=400.0,\
            dNdE_option="1"):

        if dNdE_option=="1":
            self.integrate_dNdE=integrate_dNdE_1
        else:
            self.integrate_dNdE=integrate_dNdE_1
        self.Mchi=Mchi #GeV #TODO: fill in realistic value for mass

        #TODO: fill in realistic values for cross section
        #self.Slist=[0,1,2,3,4,5,6,7,8,9,10]
        #self.Slist=[10**(-26+(26-19)*n/100) for n in range(101)]
        self.Slist=constants.Slist

        self.observations=[]
        self.P_S_sum=0
        self.memoizedict={}
        self.memoizedict2={}

    def load_observations(self, filelist):
        f=open(filelist, 'rb')
        self.observations=[]
        for line in f.readlines():
            filename=line.strip()
            o=Observation()
            o=o.open(filename)
            self.observations.append(o)

    def P_Nbpi(self, obs, Nbpi):
        top=PoissonApprox2(obs.Nmi,obs.Ci*Nbpi)

        if obs.Nmi>5:
            bottom=1.0/(obs.Ci)
        else:
            #memoize the bottom part
            if obs.Nmi in self.memoizedict2:
                if obs.Ci*Nbpi in self.memoizedict2[obs.Nmi]:
                    bottom=self.memoizedict2[obs.Nmi][obs.Ci*Nbpi]
                else:
                    self.memoizedict2[obs.Nmi][obs.Ci*Nbpi]=math.fsum(\
                        [PoissonApprox2(obs.Nmi,obs.Ci*Nbpi2)\
                        for Nbpi2 in range(obs.Nmi*3)])
                    bottom=self.memoizedict2[obs.Nmi][obs.Ci*Nbpi]
            else:
                self.memoizedict2[obs.Nmi]={}
                self.memoizedict2[obs.Nmi][obs.Ci*Nbpi]=math.fsum(\
                    [PoissonApprox2(obs.Nmi,obs.Ci*Nbpi2)\
                    for Nbpi2 in range(obs.Nmi*3)])
                bottom=self.memoizedict2[obs.Nmi][obs.Ci*Nbpi]
        prob=top/bottom

        return prob

    def P_Npi_JiAiEiNbpiS(self, obs, iJi, iAi, iEi, Nbpi, S):
        Mchi=self.Mchi
        Ngi=(S/(8*math.pi*(Mchi**2)))*\
            self.integrate_dNdE(obs.Eirange[iEi][0],obs.Eirange[iEi][1])\
            *obs.Ai[iEi][iAi]*obs.Ti*obs.Ji[iJi]
        lmbda=Nbpi+Ngi
        k=obs.Npi
        prob=PoissonApprox2(k,lmbda)
        #print("P_Npi_JiAiEiNbpiS: %s" % prob)
        return PoissonApprox2(k,lmbda)

    def P_Npi_NbpiS(self, obs, Nbpi, S):
        p=[]
        for iEi, Ei in enumerate(obs.Ei):
            for iAi, Ai in enumerate(obs.Ai[iEi]):
                for iJi, Ji in enumerate(obs.Ji):
                    P=self.P_Npi_JiAiEiNbpiS(obs, iJi, iAi, iEi, Nbpi, S)
                    P=P*obs.P_Ji[iJi]*obs.P_Ai[iEi][iAi]*obs.P_Ei[iEi]
                    p.append(P)
        prob=math.fsum(p)
        #print("P_Npi_NbpiS: %s" % prob)
        return math.fsum(p)

    def P_Npi_S(self, obs, S):
        Nbpi_cutoff=int(3*obs.Nmi/obs.Ci)
        if Nbpi_cutoff<30:
            Nbpi_cutoff=30
        p=[self.P_Nbpi(obs, Nbpi)*self.P_Npi_NbpiS(obs, Nbpi, S)\
            for Nbpi in range(Nbpi_cutoff)]
        prob=math.fsum(p)
        #print("P_Npi_S: %s" % prob)
        return prob

    def P_S_Np_product(self, S):
        Mchi=self.Mchi
        T=[]
        for obs in self.observations:
            P=self.P_Npi_S(obs, S)
            if P!=0.0:
                T.append(math.log(P))
            else:
                return 0.0
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
    
    def P_S_Np_all(self):
        P_Slist=[]
        for j,S in enumerate(self.Slist):
            P_Slist.append(self.P_S_Np(self.Slist[j]))
        return P_Slist

def integrate_dNdE_1(Ei0, Ei1):
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
    parser = argparse.ArgumentParser(description='filelist is data, Mchi is in GeV, and dNdE_option is which channel')
    parser.add_argument('--filelist', type=str, nargs=1, required=True)
    parser.add_argument('--Mchi', type=float, nargs=1, required=True)
    parser.add_argument('--dNdE_option', type=str, nargs=1, required=True)
    args=parser.parse_args()
    Mchi=args.Mchi[0]
    dNdE_option=args.dNdE_option[0]
    filelist=args.filelist[0]
    print("Options: filelist=%s, Mchi=%s, dNdE_option=%s"\
        % (filelist, Mchi, dNdE_option))
    print("Loading data...")
    b=Bayes(Mchi=Mchi, dNdE_option=dNdE_option)
    b.load_observations(filelist)
    print("Calculating probabilities...")
    Slist=b.Slist
    P_Slist=b.P_S_Np_all()
    P_sum=0
    print("S\tP\tP_sum")
    for j,S in enumerate(Slist):
        P=P_Slist[j]
        P_sum=P_sum+P
        print("%s\t%s\t%s" % (S, P, P_sum))
    #print("Sum over P: %s" % math.fsum(\
    #   [b.P_S_Np(S) for S in b.Slist]))
