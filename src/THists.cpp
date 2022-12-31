#include "THists.h"
#include "TH3.h"
#include "TMath.h"
#include "TProfile.h"
#include "TString.h"

THists::THists(std::string& name)
{
    histfile = TFile::Open(name.c_str(), "recreate");
    float a[11]{0.16, 0.316, 0.631, 1.26, 2.5, 5, 10, 20, 40, 80 ,160};
    for(int i = 0; i < 10; i++)
    {
        h100[i] = new TH2D(Form("h%d", i+100), Form("Core position of Event Triggered >= 4tel .%f - %f TeV", a[i], a[i+1]), 300, -1500, 1500, 300, -1500, 1500);
        h200[i] = new TProfile(Form("h%d", i+200), Form("Image Dist Versus Rp (%f - %f TeV)", a[i], a[i+1]), 50, 0, 1000, 0, 5);
        h250[i] = new TProfile(Form("h%d", i + 250), Form("Beta versus Dist (%f - %f TeV)", a[i],a[i+1]), 500, 0, 5,  0, 5);
        h2250[i] = new TProfile( Form("h%d ", i+ 2250), Form("MISS Versus Dist (%f - %f TeV)", a[i], a[i + 1]),  500, 0, 5,  0, 5);
        h400[i] = new TH1D(Form("h%d", i+400), Form("MRSW Distribution (%f - %f TeV)", a[i], a[i+1]), 400, -10, 10);
        h450[i] = new TProfile(Form("h%d", i+450), Form("Width Distribution Versus Rp (%f - %f TeV)", a[i], a[i+1]), 50, 0, 1000, 0, 4);
        h550[i] = new TProfile(Form("h%d", i+550), Form("Length Distribution Versus Rp (%f - %f TeV)", a[i], a[i+1]), 50, 0, 1000,  0,4);
        h500[i] = new TH1D(Form("h%d", i+500), Form("MRSL Distribution (%f - %f TeV)", a[i], a[i+1]), 400, -10, 10);
        h600[i] = new TProfile(Form("h%d", i+600), Form("Direction Error versus Core array distance (%f - %f TeV)", a[i], a[i+1]), 50, 0, 1000,  0, 10);
    }
        h100[10] = new TH2D(Form("h%d", 10+100), Form("Core position of Event Triggered >= 4tel . >%f TeV", a[10]), 300, -1500, 1500, 300, -1500, 1500);
        h200[10] = new TProfile(Form("h%d", 10+200), Form("image Dist Versus Rp ( >%f  TeV)", a[10]), 50, 0, 1000,  0, 8);
        h250[10] = new TProfile(Form("h%d", 10 + 250), Form("Beta versus Dist ( >%fTeV)", a[10]), 500, 0, 5,  0, 5);
        h2250[10] = new TProfile(Form("h%d", 10 + 2250), Form("MISS versus Dist ( >%fTeV)", a[10]), 500, 0, 5, 0, 5);
        h450[10] = new TProfile(Form("h%d", 10+410), Form("Width Distribution Versus Rp ( >%f  TeV)", a[10]), 50, 0, 1000,  0, 4);
        h550[10] = new TProfile(Form("h%d", 10+510), Form("Length Distribution Versus Rp ( >%f  TeV)", a[10]),50, 0, 1000, 0, 4);
        h400[10] = new TH1D(Form("h%d", 10+400), Form("MRSW Distribution (> %f TeV)", a[10]), 200, -10, 10);
        h500[10] = new TH1D(Form("h%d", 10+500), Form("MRSL Distribution (> %f TeV)", a[10]), 200, -10, 10);
        h600[10] = new TProfile(Form("h%d", 10+600), Form("Direction Error versus Core array distance (> %f TeV)", a[10]), 50, 0, 1000,  0, 10);

        h1 = new TH1D("h1", "True Energy Histogram (All events)", 20, -1, 3);
        h2 = new TH1D("h2", "True Energy Histogram (Trigger events)", 20, -1, 3);
        h3 = new TH1D("h3", "True Energy Histogram (>4 Image events)", 20, -1, 3);
        h4 = new TH2D("h4", "Core Array Distance versus True Energy (>4 Image events)", 100, 0, 1000, 20, -1, 3);
        h5 = new TH1D("h5", "True Energy Histogram(Pass Shape Cut)" , 20, -1, 3);
        h6 = new TH1D("h6", "True Energy Histogram(Pass theta2 Cut)" , 20, -1, 3);
        h301 = new TH2D("h301", "MRSL:MRSW ", 200, -10, 10, 200, -10, 10 );
        h302 = new TH2D("h302", "MRSL:MRSW direction_error < 1deg ", 200, -10, 10, 200, -10, 10 );
        h303 = new TH1D("h303", "MRSW distribution", 200, -10, 10 );
        h304 = new TH1D("h304", "MRSL distribution", 200, -10, 10 );
        h60 = new TH1D("h60", "Direction Error distribution", 2000, 0, 10 );
        h61 = new TH1D("h61", "Direction Error distribution(Pass Shape Cut)", 2000, 0, 10 );
        h62 = new TProfile("h62", "direction_error versus Energy", 20, -1, 3,  0, 2);
        h63 = new TProfile("h63", "direction_error versus Energy (Pass Shpae Cuts)", 20, -1, 3,  0, 100);
        h64 = new TH1D("h64", "theta square distribution ", 1000, 0, 0.1);
        hist_max = new TH3D("hist_max", "x:Xmax y:Rp z:Dist ",400, 0, 1200, 70, 0, 700, 100, 0, 5);
}

THists::~THists()
{
    //delete []h100;
    //delete []h200;
    //delete []h400;
    //delete []h500;
    //delete []h600;
    //delete h1, h2, h3 ,h4;
    //delete h301, h302, h303, h304;
    //delete h60, h61;

}

void THists::Write()
{
    histfile->cd();
    for(int i = 0; i < 11; i++)
    {
        h100[i]->Write();
        h200[i]->Write();
        h250[i]->Write();
        h2250[i]->Write();
        h400[i]->Write();
        h450[i]->Write();
        h500[i]->Write();
        h550[i]->Write();
        h600[i]->Write();
    }
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    h5->Write();
    h6->Write();
    h301->Write();
    h302->Write();
    h303->Write();
    h304->Write();
    h60->Write();
    h61->Write();
    h62->Write();
    h63->Write();
    h64->Write();
    hist_max->Write();
    histfile->Write();
    histfile->Close();

}