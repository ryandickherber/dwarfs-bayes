from ROOT import TCanvas, TH1F, TSlider
from ROOT import TFile, TTree

#char * filename="/data/rcdickherber/SEGUE2/61533_config_Eall/results_s6.root";
filename="/data/rcdickherber/SEGUE2/61533_config_Eall/results_s6.root"
#TFile * tf=new TFile(filename);
tf=TFile(filename)
#TTree * eventTree=(TTree*)tf->Get("EventStatsTree");
eventTree=tf.Get("EventStatsTree")

#UInt_t ArrayEventNum;
#eventTree->SetBranchAddress("ArrayEventNum", &ArrayEventNum);
#UInt_t RunNum;
#eventTree->SetBranchAddress("RunNum", &RunNum);
#UInt_t MJDInt;
#eventTree->SetBranchAddress("MJDInt", &MJDInt);
#Double_t MJDDbl;
#eventTree->SetBranchAddress("MJDDbl", &MJDDbl);
#Bool_t OnEvent;
#eventTree->SetBranchAddress("OnEvent", &OnEvent);
#Bool_t OffEvent;
#eventTree->SetBranchAddress("OffEvent", &OffEvent);
#Double_t Azimuth;
#eventTree->SetBranchAddress("Azimuth", &Azimuth);
#Float_t Elevation;
#eventTree->SetBranchAddress("Elevation", &Elevation);
#Float_t EnergyGeV;
#eventTree->SetBranchAddress("EnergyGeV", &EnergyGeV);
#Float_t EnergyRMS;
#eventTree->SetBranchAddress("EnergyRMS", &EnergyRMS);
#Float_t Noise;
#eventTree->SetBranchAddress("Noise", &Noise);
#Float_t Offset;
#eventTree->SetBranchAddress("Offset", &Offset);
#Float_t EffectiveArea;
#eventTree->SetBranchAddress("EffectiveArea", &EffectiveArea);
#
#
#int entries=eventTree->GetEntries();

for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
    print eventTree.ArrayEventNum

#
#for (int i=0; i<entries; i++) {
#    eventTree->GetEntry(i);
#    if (OnEvent || OffEvent) {
#        cout << RunNum << " " << ArrayEventNum << " " << EffectiveArea << endl;
#   //if (EffectiveArea>10e-3)
#   //    cout << EffectiveArea << endl;
#    }
#}
