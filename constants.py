import math

NE=31
Estart=100.0
Eincrement=0.15
Erange=[]
for i in range(NE):
    if (i==0):
        Erange.append(Estart)
    else:
        Erange.append(Erange[i-1]+Eincrement*Erange[i-1])

NM=30
Mstart=115.0
Mincrement=0.15
Mlist=[]
for i in range(NM):
    if (i==0):
        Mlist.append(Mstart)
    else:
        Mlist.append(Mlist[i-1]+Mlist[i-1]*Mincrement)


#HACK alert! Sadly, this number is necessarily.
#It is canceled, but if I had to write this software again it wouldn't
#be necessary. It becomes necessary to set this number arbitrarily high
#when calculating lots of data due to the smallness of the P_Npi_S
#probabilities.
P_Npi_S_factor=100

#Slist=[0,1,2,3,4,5,6,7,8,9,10]
#Slist=[10**(-26+(26-19)*n/100) for n in range(101)]
#Slist=[math.pow(10,(-29+(29-21)*(1.0*n)/100)) for n in range(101)]
#Slist=[math.pow(10,(-30+(30-17)*(1.0*n)/100)) for n in range(101)]

#use this:
#Slist=[math.pow(10,(-40+(40-19)*(1.0*n)/100)) for n in range(101)]
NS=100
f=lambda n : math.pow(10,(-26+(26-21)*(1.0*n)/NS))
Slist=[f(n) for n in range(NS+1)]
Slist_sum=math.fsum(Slist)

#or...
#start=10e-26
#finish=10e-22
#increment=(finish-start)/100
#Slist=[start+i*increment for i in range(100)]

#TODO: Is this right?
Slist_diff=[f(n+1)-f(n) for n in range(NS+1)]
Slist_integration_factors={}
for n in range(NS+1):
    Slist_integration_factors[f(n)]=Slist_diff[n]/math.fsum(Slist_diff)

#J-factors... units are GeV^2cm^-5
J_SEGUE=7.7e17
J_Draco=4*3.832e17
J_UrsaMinor=7*3.832e17
J_BOOTES1=3*3.832e17
J_WilmanI=22*3.832e17

