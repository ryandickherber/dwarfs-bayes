NE=31
Estart=200.0
Eincrement=0.15
Erange=[]
for i in range(NE):
    if (i==0):
        Erange.append(Estart)
    else:
        Erange.append(Erange[i-1]+Eincrement*Erange[i-1])

NM=30
Mstart=200.0
Mincrement=0.15
Mlist=[]
for i in range(NM):
    if (i==0):
        Mlist.append(Mstart)
    else:
        Mlist.append(Mstart+Mlist[i-1]*Mincrement)

#Slist=[0,1,2,3,4,5,6,7,8,9,10]
Slist=[10**(-26+(26-19)*n/100) for n in range(101)]

#J-factors... units are GeV^2cm^-5
J_SEGUE=7.7e17
J_Draco=4*3.832e17
J_UrsaMinor=7*3.832e17
J_BOOTES1=3*3.832e17
J_WilmanI=22*3.832e17
