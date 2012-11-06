from ROOT import TCanvas, TH1F, TSlider
from ROOT import TFile, TTree
import bayes
import constants
import math

#observations=[]

#ArrayEventNum
#RunNum
#MJDInt
#MJDDbl
#OnEvent
#OffEvent
#Azimuth
#Elevation
#EnergyGeV
#EnergyRMS
#Noise
#Offset
#EffectiveArea

myf=open("data_test_5.stage6")
filenames=[x.strip() for x in myf.readlines()]
myf.close()

sourcelist=["SEGUE","Draco","UrsaMinor","WilmanI","BOOTES1"]

for source in sourcelist:
    for i in range(constants.NE-1):
        Erange0=constants.Erange[i]
        Erange1=constants.Erange[i+1]
        #First, get all alpha...
        alphalist=[]
        for filename in filenames:
            if source not in filename:
                continue
            tf=TFile(filename)
            runTree=tf.Get("RunStatsTree")
            runTree.GetEntry(0)
            alpha=runTree.fAlpha
            if alpha not in alphalist:
                alphalist.append(alpha)
        for alpha in alphalist:
            obs=bayes.Observation()
            obs.Ci=1/alpha
            obs.Ti=0
            obs.EnergyBin=[Erange0,Erange1]
            obs.EnergyBinIndex=i
            obs.EnergyBinIndexes=[i]
            obs.Npi=0
            obs.Nmi=0
            areas=[]
            energies=[]
            for filename in filenames:
                if source not in filename:
                    continue
                tf=TFile(filename)
                runTree=tf.Get("RunStatsTree")
                runTree.GetEntry(0)
                if alpha!=runTree.fAlpha:
                    #print("%s!=%s" % (alpha, runTree.fAlpha))
                    continue
                obs.Ti=obs.Ti+runTree.faLiveTime
                eventTree=tf.Get("EventStatsTree")
                for j in range(eventTree.GetEntries()):
                    eventTree.GetEntry(j)
                    if eventTree.EnergyGeV>=Erange0\
                        and eventTree.EnergyGeV<Erange1\
                        and (eventTree.OnEvent or eventTree.OffEvent):
                        if eventTree.OnEvent:
                            obs.Npi=obs.Npi+1
                        if eventTree.OffEvent:
                            obs.Nmi=obs.Nmi+1
                        areas.append(eventTree.EffectiveArea)
                        energies.append(eventTree.EnergyGeV)
            nevents=len(energies)
            print("Ci: %s, EnergyBinIndex: %s, nevents: %s" %
                (obs.Ci, i, nevents))
            if nevents<10:
                continue
            avgareas=math.fsum(areas)/len(areas)
            obs.Ai=[[avgareas]]
            obs.P_Ai=[[1.0]]
            avgenergies=math.fsum(energies)/len(energies)
            obs.Ei=[avgenergies]
            obs.Eirange=[(Erange0,Erange1)]
            obs.P_Ei=[1.0]
            obs.source=source
            if source=="SEGUE":
                obs.Ji=[constants.J_SEGUE]
            elif source=="Draco":
                obs.Ji=[constants.J_Draco]
            elif source=="UrsaMinor":
                obs.Ji=[constants.J_UrsaMinor]
            elif source=="WilmanI":
                obs.Ji=[constants.J_WilmanI]
            elif source=="BOOTES1":
                obs.Ji=[constants.J_BOOTES1]
            else:
                print("ERROR! No J-factor found...")
                exit()
            obs.P_Ji=[1.0]
            obs.save("data_test_5/%s_%s_%s.pickle"\
                % (obs.EnergyBinIndex,obs.Ci,obs.source))
