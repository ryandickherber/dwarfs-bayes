
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
        Mlist.append(Mstart+Mlist[i]*Mincrement)
