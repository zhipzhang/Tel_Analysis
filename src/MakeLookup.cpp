/*
        Generate the Lookup Table For Reconstruction 
        Width / Length
        Energy

*/
#include "TLookup_table.h"
#include <string>
#include "TImage_Parameter.h"
#include "TRecData.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"



void syntax();
void syntx()
{
    printf("Usage: MakeLookup output_file input_file \n");
}


int main(int argc, char** argv)
{
    const char* input_file = NULL;
    std::string out_file = "lookup_table.root";
    TImage_Parameter* image = new TImage_Parameter();
    TRecData* rec = new TRecData();
    TLookup_table* lookup = new TLookup_table();
    
    while(argc > 2)
    {
        if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") ==0)
        {
            syntax();
            argc --;
            argv ++;
            continue;
        }
        if(strcmp(argv[1], "--out_file") == 0)
        {
            out_file  = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        else
        {
            break;
        }

    }
    TFile* out_root = TFile::Open(out_file.c_str(),"recreate");
    TH2D * h1 = new TH2D("h1","Rp versus log10A (no weighted)",100, 0, 450, 80, 0, 8);
    TH2D * h10 = new TH2D("h10","Rp versus log10A (weighted)",100, 0, 450, 80, 0, 8);
    TH2D * h2 = new TH2D("h2","Rp versus log10A (weighted by E/A)",100, 0, 450, 80, 0, 8);
    TH2D * h22 = new TH2D("h22","Rp versus log10A (weighted by square(E/A))",100, 0, 450, 80, 0, 8);
    TH2D * h3 = new TH2D("h3","Rp versus log10A (weighted by length)",100, 0, 450, 80, 0, 15);
    TH2D * h33 = new TH2D("h33","Rp versus log10A (weighted by  square length)",100, 0, 450, 80, 0, 15);
    TH2D * h4 = new TH2D("h4","Rp versus log10A (weighted by width)",100, 0, 450, 80, 0, 15);
    TH2D * h44 = new TH2D("h44","Rp versus log10A (weighted by square width)",100, 0, 450, 80, 0, 15);
    while(argc > 1 || input_file != NULL )
    {
        input_file = argv[1];
        argc--;
        argv++;
        TFile* input_root = TFile::Open(input_file, "read");
        input_file = NULL;
        TTree* image_tree = (TTree*) input_root->Get("image_parameter");
        TTree* rec_tree = (TTree*) input_root->Get("rec_data");

        image_tree->SetBranchAddress("image_parameter", &image);
        rec->InitRead(rec_tree);
        for(int i = 0 ; i< image_tree->GetEntries(); i++)
        {
            image_tree->GetEntry(i);
            rec_tree->GetEntry(i);
            for(int j = 0; j < image->image_tel.size(); j++)
            {
                int itel = image->image_tel[j];
                double size = image->GetTelSize(itel);
                double rp = image->GetTelRp(itel);
                double weight = rec->weight;
                if(size > 100 && rp < 400)
                {
                    h1->Fill(log10(size), rp);
                    h10->Fill(log10(size), rp, weight);
                    h2->Fill(log10(size), rp, rec->GetEnergy()/size * weight);
                    h22->Fill(log10(size), rp, pow(rec->GetEnergy()/size, 2) * weight);
                    h3->Fill(log10(size), rp, image->GetTelLength(itel) * weight);
                    h33->Fill(log10(size), rp, pow(image->GetTelLength(itel), 2) * weight);
                    h4->Fill(log10(size), rp, image->GetTelwidth(itel) * weight);
                    h44->Fill(log10(size), rp, pow(image->GetTelwidth(itel), 2) * weight);
                }

            }
        }
        out_root->cd();

        TH2D* he = new TH2D("he", "mean_energy_lookup", 100, 0, 450, 80, 0, 8);
        TH2D* hee = new TH2D("hee", "var_energy_lookup", 100, 0, 450, 80, 0, 8);
        TH2D* hw = new TH2D("hw", "mean_width_lookup", 100, 0, 450, 80, 0, 8);
        TH2D* hww = new TH2D("hww", "var_width_lookup", 100, 0, 450, 80, 0, 8);
        TH2D* hl = new TH2D("hl", "mean_length_lookup", 100, 0, 450, 80, 0, 8);
        TH2D* hll = new TH2D("hll", "var_length_lookup", 100, 0, 450, 80, 0, 8);

        for(int i = 0; i < h1->GetXaxis()->GetNbins(); i++)
        {
            for(int j = 0; j < h1->GetYaxis()->GetNbins(); j++)
            {
                he->Fill(h1->GetXaxis()->GetBinCenter(i), h1->GetYaxis()->GetBinCenter(j), h2->GetBinContent(i,j)/h10->GetBinContent(i, j));
                hee->Fill(h1->GetXaxis()->GetBinCenter(i), h1->GetYaxis()->GetBinCenter(j), h22->GetBinContent(i,j)/h10->GetBinContent(i, j));
                hw->Fill(h1->GetXaxis()->GetBinCenter(i), h1->GetYaxis()->GetBinCenter(j), h3->GetBinContent(i,j)/h10->GetBinContent(i, j));
                hww->Fill(h1->GetXaxis()->GetBinCenter(i), h1->GetYaxis()->GetBinCenter(j), h33->GetBinContent(i,j)/h10->GetBinContent(i, j));
                hl->Fill(h1->GetXaxis()->GetBinCenter(i), h1->GetYaxis()->GetBinCenter(j), h4->GetBinContent(i,j)/h10->GetBinContent(i, j));
                hll->Fill(h1->GetXaxis()->GetBinCenter(i), h1->GetYaxis()->GetBinCenter(j), h44->GetBinContent(i,j)/h10->GetBinContent(i, j));

            }
        }

        out_root->cd();
        he->Write();
        hee->Write();
        hw->Write();
        hww->Write();
        hl->Write();
        hll->Write();



    }

    


    

}