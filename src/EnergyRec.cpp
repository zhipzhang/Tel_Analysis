/*
        
    A small program to rec energy, In Fact main program Analysis can Also do the Energy Rec.
    This Program can Only Generate a histogram of E_true and E_rec


*/
#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"

#include "TH2D.h"
#include "TRecData.h"
#include "TImage_Parameter.h"



int main(int argc, char** argv)
{
    std::string outfile = "rec_energy.root";
    std::string lookupfile = "";
    std::string inputfile = "";
    TImage_Parameter* image = new TImage_Parameter();
    TRecData* rec = new TRecData();

    while(argc > 2)
    {
        if(strcmp(argv[1], "--out-file") == 0)
        {
            outfile = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        if(strcmp(argv[1], "--lookup-file") == 0)
        {
            lookupfile = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        else
        {
            break;
        }
    }

    if(lookupfile.empty())
    {
        std::cout << "No Lookup File! Program Exit!" << std::endl;
        exit(EXIT_FAILURE);
    }
    TFile* out_root = TFile::Open(outfile.c_str(), "recreate");
    TH2D* h1 = new TH2D("h1", "Rec Energy Verus True Energy", 30, 0, 3, 600, -1, 3);
    TFile* lookup_root = TFile::Open(lookupfile.c_str(), "read");
    if(lookup_root->IsZombie())
    {
        std::cout << "Can't Open Lookup File !" << std::endl;
        exit(EXIT_FAILURE);
    }
    TH2D *he = (TH2D*) ( lookup_root->Get("he"));
    TH2D* hee = (TH2D*)lookup_root->Get("hee");
    TH2D* hl = (TH2D*)lookup_root->Get("hl");
    TH2D* hw = (TH2D*)lookup_root->Get("hw");
    TH2D* hll = (TH2D*)lookup_root->Get("hll");
    TH2D* hww = (TH2D*)lookup_root->Get("hww");
    
    while(argc > 1)
    {
        inputfile = argv[1];
        argc --;
        TFile* input_root = TFile::Open(inputfile.c_str(), "read");
        TTree* rec_tree = (TTree*)input_root->Get("rec_data");
        TTree* image_tree = (TTree*) input_root->Get("image_parameter");
        image_tree->SetBranchAddress("image_parameter", &image);
        rec->InitRead(rec_tree);

        for(int i = 0; i < image_tree->GetEntries(); i++)
        {
            image_tree->GetEntry(i);
            rec_tree->GetEntry(i);
            double all_e = 0;
            double all_w = 0;
            for(int j = 0; j < image->image_tel.size(); j++)
            {
                int itel = image->image_tel[j];
                double size = image->GetTelSize(itel);
                double rp = image->GetTelRp(itel);
                if(size < 100)
                {
                    continue;
                }
                int xbin = he->GetXaxis()->FindBin(rp);
                int ybin = he->GetYaxis()->FindBin(log10(size));
                double energy = size/he->GetBinContent(xbin, ybin) ;
                double var_e;
                if(he->GetBinContent(xbin, ybin) > 0)
                    var_e =  hee->GetBinContent(xbin, ybin)/he->GetBinContent(xbin, ybin);
                else
                var_e = 999.;
                double weight = 1/(0.01 + var_e * var_e);
                all_e += log(energy) * weight;
                all_w += weight;
             }
             if(all_w != 0)
             {
                double rec_energy = exp(all_e/all_w);
                h1->Fill(rec->GetEnergy(), log10(exp(all_e/all_w)/rec->GetEnergy()), pow(rec->GetEnergy(), -1));
             }
             else
             {
                double ggg = rec->GetEnergy();
             }

        }
        input_root->Close();
    } 
    out_root->cd();
    h1->Write();
    out_root->Close();


    }


