from ROOT import TCanvas, TH1F, TSlider
from ROOT import TFile, TTree
import bayes
import constants
import math

observations=[]

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

filenames=[\
"/data/rcdickherber/SEGUE2/61533_config_Eall/results_s6.root",\
]
#for filename in filenames:
#    tf=TFile(filename)
#    eventTree=tf.Get("EventStatsTree")
#    eventTree.GetEntry(0)
#    print(eventTree.RunNum)

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

            #for the time being, drop observations that have no on or off events
            if obs.Npi==0 and obs.Nmi==0:
                continue

            avgareas=math.fsum(areas)/len(areas)
            obs.Ai=[[avgareas]]
            obs.P_Ai=[[1.0]]
            avgenergies=math.fsum(energies)/len(energies)
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
            obs.Ci=runstatsTree.faLiveTime
            obs.Ti=runstatsTree.fAlpha
            observations.append(obs)

for obs in observations:
    #print(obs.Npi, obs.Nmi)
    print obs

#
#for (int i=0; i<entries; i++) {
#    eventTree->GetEntry(i);
#    if (OnEvent || OffEvent) {
#        cout << RunNum << " " << ArrayEventNum << " " << EffectiveArea << endl;
#   //if (EffectiveArea>10e-3)
#   //    cout << EffectiveArea << endl;
#    }
#}
