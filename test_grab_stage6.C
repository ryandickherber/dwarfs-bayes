#include <vector>
#include <string>
#include <cstring>
#include <iostream>
using namespace std;


int test_grab_stage6()
{

char * filename="/data/rcdickherber/SEGUE2/61533_config_Eall/results_s6.root";
TFile * tf=new TFile(filename);
TTree * eventTree=(TTree*)tf->Get("EventStatsTree");

UInt_t ArrayEventNum;
eventTree->SetBranchAddress("ArrayEventNum", &ArrayEventNum);
UInt_t RunNum;
eventTree->SetBranchAddress("RunNum", &RunNum);
UInt_t MJDInt;
eventTree->SetBranchAddress("MJDInt", &MJDInt);
Double_t MJDDbl;
eventTree->SetBranchAddress("MJDDbl", &MJDDbl);
Bool_t OnEvent;
eventTree->SetBranchAddress("OnEvent", &OnEvent);
Bool_t OffEvent;
eventTree->SetBranchAddress("OffEvent", &OffEvent);
Double_t Azimuth;
eventTree->SetBranchAddress("Azimuth", &Azimuth);
Float_t Elevation;
eventTree->SetBranchAddress("Elevation", &Elevation);
Float_t EnergyGeV;
eventTree->SetBranchAddress("EnergyGeV", &EnergyGeV);
Float_t EnergyRMS;
eventTree->SetBranchAddress("EnergyRMS", &EnergyRMS);
Float_t Noise;
eventTree->SetBranchAddress("Noise", &Noise);
Float_t Offset;
eventTree->SetBranchAddress("Offset", &Offset);
Float_t EffectiveArea;
eventTree->SetBranchAddress("EffectiveArea", &EffectiveArea);


int entries=eventTree->GetEntries();

for (int i=0; i<entries; i++) {
    eventTree->GetEntry(i);
    if (OnEvent || OffEvent) {
        cout << RunNum << " " << ArrayEventNum << " " << EffectiveArea << endl;
	//if (EffectiveArea>10e-3)
	//    cout << EffectiveArea << endl;
    }
}

return 0;
}

//TTree * tt2=(TTree*)tf->Get("ParameterisedEvents/ParEventsTree");
//VAShowerData * sd=new VAShowerData();
//VAHillasData * hd=new VAHillasData();
//tt1->SetBranchAddress("S",&sd);
//tt2->SetBranchAddress("H",&hd);
//
//int entries=tt1->GetEntries();
//for (int entry=0; entry<entries; entry++)
//{
//	//whole array parameters
//	tt1->GetEntry(entry);
//	double fMSW=sd->fMSW;
//	double fMSL=sd->fMSL;
//	double fArrayEventNum=sd->fArrayEventNum;
//	double fEnergy_GeV=sd->fEnergy_GeV;
//
//	//single telescope parameters
//	//assuming exactly 4 scopes per full-array entry
//	int hfTelId[4];
//	int hfArrayEventNum[4];
//	int hfPixelsInImage[4];
//	for (int myit=0; myit<4; myit++)
//	{
//		tt2->GetEntry(4*entry+myit);
//		int root_fucking_sucks=hd->fTelId;
//		hfTelId[myit]=root_fucking_sucks;
//		hfArrayEventNum[myit]=hd->fArrayEventNum;
//		hfPixelsInImage[myit]=hd->fPixelsInImage;
//	}
//
//	/*
//	if (!(hfTelId[0]==0 && hfTelId[1]==1 && hfTelId[2]==2 && hfTelId[3]==3))
//	{
//		cout << "Fatal error! Telescopes don't match up" << endl;
//		exit();
//	}
//	*/
//
//	int energy_bin;
//	//cout << "fEnergyGeV: " << fEnergy_GeV << endl;
//	energy_bin=get_energy_bin(fEnergy_GeV); //reconstructed energy
//
//	/*
//
//	//quality cut
//	int cut_hfPixelsInImage=3;
//	if (hfPixelsInImage[0]<cut_hfPixelsInImage
//		|| hfPixelsInImage[1]<cut_hfPixelsInImage
//		|| hfPixelsInImage[2]<cut_hfPixelsInImage
//		|| hfPixelsInImage[3]<cut_hfPixelsInImage)
//	{
//		//cout << "quality cut" << endl;
//		continue;
//	}
//	else
//	{
//		//cout << "quality pass" << endl;
//	}
//	*/
//
//	double electron_cut_param=1.0*fMSW/4+0.0*fMSL/4;
//
//	//electron-like cut
//	/*
//	if (electron_cut_param<1.3/4)
//	{
//		//cout << "electron cut" << endl;
//		continue;
//	}
//	else
//	{
//		//cout << "electron pass" << endl;
//	}
//	//cout << "both cuts pass" << endl;
//	//cout << "msw: " << fMSW << endl;
//	*/
//
//	hist[energy_bin]->Fill(electron_cut_param);
//}
//
//delete sd;
//delete hd;
//delete tf;
