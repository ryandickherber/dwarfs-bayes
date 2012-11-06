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

myf=open("data_test_3.stage6")
filenames=[x.strip() for x in myf.readlines()]
myf.close()

for filename in filenames:
    tf=TFile(filename)
    eventTree=tf.Get("EventStatsTree")
    for i,e in enumerate(constants.Erange):
        if i<constants.NE-1:
            Erange0=constants.Erange[i]
            Erange1=constants.Erange[i+1]
            eventTree.GetEntry(0)
            obs=bayes.Observation()
            obs.RunNum=eventTree.RunNum
            obs.EnergyBin=[Erange0, Erange1]
            obs.EnergyBinIndex=i
            obs.EnergyBinIndexes=[i]
            obs.Npi=0
            obs.Nmi=0
            areas=[]
            energies=[]
            for i in range(eventTree.GetEntries()):
                eventTree.GetEntry(i)
                if eventTree.EnergyGeV>=Erange0\
                    and eventTree.EnergyGeV<=Erange1\
                    and (eventTree.OnEvent or eventTree.OffEvent):
                    if eventTree.OnEvent:
                        obs.Npi=obs.Npi+1
                    if eventTree.OffEvent:
                        obs.Nmi=obs.Nmi+1
                    areas.append(eventTree.EffectiveArea)

                    #for now we treat reconstructed energy as real energy
                    energies.append(eventTree.EnergyGeV)
            if len(areas)!=0:
                avgareas=math.fsum(areas)/len(areas)
            else:
                #TODO: Fix this.
                avgareas=1.0
            obs.Ai=[[avgareas]]
            obs.P_Ai=[[1.0]]
            if len(energies)!=0:
                avgenergies=math.fsum(energies)/len(energies)
            else:
                #TODO: Fix this.
                avgenergies=1.0
            obs.Ei=[avgenergies]
            obs.Eirange=[(Erange0,Erange1)]
            obs.P_Ei=[1.0]
            if "SEGUE" in filename:
                obs.Ji=[constants.J_SEGUE]
            elif "Draco" in filename:
                obs.Ji=[constants.J_Draco]
            elif "UrsaMinor" in filename:
                obs.Ji=[constants.J_UrsaMinor]
            elif "WilmanI" in filename:
                obs.Ji=[constants.J_WilmanI]
            elif "BOOTES1" in filename:
                obs.Ji=[constants.J_BOOTES1]
            else:
                print("ERROR! No J-factor found...")
                exit()
            obs.P_Ji=[1.0]
            runstatsTree=tf.Get("RunStatsTree")
            runstatsTree.GetEntry(0)
            #obs.Ci=runstatsTree.faLiveTime
            #obs.Ti=runstatsTree.fAlpha
            obs.Ci=1/runstatsTree.fAlpha
            obs.Ti=runstatsTree.faLiveTime
            #observations.append(obs)
            print("%s\t%s" % (obs.RunNum,obs.EnergyBinIndex))
            obs.save("data_test_3/%s_%s.pickle"\
                % (obs.RunNum,obs.EnergyBinIndex))

