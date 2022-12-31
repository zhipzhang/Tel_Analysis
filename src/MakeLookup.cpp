/*
        Generate the Lookup Table For Reconstruction 
        Width / Length
        Energy

*/
//#include "TLookup_table.h"
#include <string>
#include "TImage_Parameter.h"
#include "TUserCuts.h"
#include "TRecData.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"



void syntax();
void syntax()
{
    printf("Usage: MakeLookup output_file input_file \n");
}


int main(int argc, char** argv)
{
    const char* input_file = NULL;
    std::string out_file = "lookup_table.root";
    TImage_Parameter* image = new TImage_Parameter();
    TRecData* rec = new TRecData();
    TUserCuts* tcuts = new TUserCuts();
   // TLookup_table* lookup = new TLookup_table();
    
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
    TH2D * h1 = new TH2D("h1","Rp versus log10A (no weighted)",50, 0, 500, 160, 0, 8);
    TH2D * sh1 = new TH2D("sh1","Rp versus log10A (no weighted smoothed) ",50, 0, 500, 160, 0, 8);
    TH2D * sh11 = new TH2D("sh11","Rp versus log10A (no weighted smoothed) For energy Rec ",50, 0, 500, 160, 0, 8);
    TH2D * h10 = new TH2D("h10","Rp versus log10A (weighted)",50, 0, 500, 160, 0, 8);
    TH2D * sh10 = new TH2D("sh10","Rp versus log10A (weighted smoothed)",50, 0, 500, 160, 0, 8);
    TH2D * sh101 = new TH2D("sh101","Rp versus log10A (weighted smoothed) For energy Rec",50, 0, 500, 160, 0, 8);
    TH2D * h2 = new TH2D("h2","Rp versus log10A (weighted by E)",50, 0, 500, 160, 0, 8);
    TH2D * sh2 = new TH2D("sh2","Rp versus log10A (weighted by E smoothed)",50, 0, 500, 160, 0, 8);
    TH2D * h22 = new TH2D("h22","Rp versus log10A (weighted by square(E))",50, 0, 500, 160, 0, 8);
    TH2D * sh22 = new TH2D("sh22","Rp versus log10A (weighted by square(E) smoothed)",50, 0, 500, 160, 0, 8);
    TH2D * h3 = new TH2D("h3","Rp versus log10A (weighted by length)",50, 0, 500, 160, 0, 8);
    TH2D * sh3 = new TH2D("sh3","Rp versus log10A (weighted by length smoothed) ",50, 0, 500, 160, 0, 8);
    TH2D * h33 = new TH2D("h33","Rp versus log10A (weighted by  square length)",50, 0, 500, 160, 0, 8);
    TH2D * sh33 = new TH2D("sh33","Rp versus log10A (weighted by  square length smoothed)",50, 0, 500, 160, 0, 8);
    TH2D * h4 = new TH2D("h4","Rp versus log10A (weighted by width)",50, 0, 500, 160, 0, 8);
    TH2D * sh4 = new TH2D("sh4","Rp versus log10A (weighted by width smoothed)",50, 0, 500, 160, 0, 8);
    TH2D * h44 = new TH2D("h44","Rp versus log10A (weighted by square width)",50, 0, 500, 160, 0, 8);
    TH2D * sh44 = new TH2D("sh44","Rp versus log10A (weighted by square width smoothed)",50, 0, 500, 160, 0, 8);
    TH2D* he = new TH2D("he", "mean_energy_lookup", 50, 0, 500, 160, 0, 8);
    TH2D* hee = new TH2D("hee", "var_energy_lookup", 50, 0, 500, 160, 0, 8);
    TH2D* hw = new TH2D("hw", "mean_width_lookup", 50, 0, 500, 160, 0, 8);
    TH2D* hww = new TH2D("hww", "var_width_lookup", 50, 0, 500, 160, 0, 8);
    TH2D* hl = new TH2D("hl", "mean_length_lookup", 50, 0, 500, 160, 0, 8);
    TH2D* hll = new TH2D("hll", "var_length_lookup", 50, 0, 500, 160, 0, 8);
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
            double weight = rec->weight;
            for(int j = 0; j < image->image_tel.size(); j++)
            {
                int itel = image->image_tel[j];
                double size = image->GetTelSize(itel);
                double rp = image->GetTelRp(itel);
                double dist = image->GetDist(itel) * TMath::RadToDeg();
                double weight = pow(rec->GetEnergy(), -1);
                if(size > tcuts->min_size && dist < tcuts->max_dist)
                {
                    h1->Fill(rp, log10(size));
                    h10->Fill(rp, log10(size), weight);
                    h2->Fill(rp, log10(size), size/rec->GetEnergy() * weight);
                    h22->Fill(rp, log10(size), pow(size/rec->GetEnergy(), 2) * weight);
                    h3->Fill(rp , log10(size), image->GetTelLength(itel) * weight);
                    h33->Fill(rp , log10(size), pow(image->GetTelLength(itel), 2) * weight);
                    h4->Fill(rp, log10(size),  image->GetTelwidth(itel) * weight);
                    h44->Fill(rp , log10(size),  pow(image->GetTelwidth(itel), 2) * weight);
                }
                else
                {
                    continue;
                }

            }
        }
        printf("nbins is %d",(int) h1->GetNbinsX());
        printf("nbins is %d",(int) h1->GetNbinsY());
        for( int i = 0; i < h1->GetNbinsX(); i++)
        {
            for(int j = 0 ; j < h1->GetNbinsY(); j++)
            {
                double n = h1->GetBinContent(i, j);
                double ns = h1->GetBinContent(i, j);
                double nw = h10->GetBinContent(i, j);
                double nws = h10->GetBinContent(i, j);
                double we = h2->GetBinContent(i , j);
                double wee = h22->GetBinContent(i, j);
                double wl = h3->GetBinContent(i, j);
                double wll = h33->GetBinContent(i, j);
                double ww = h4->GetBinContent(i, j);
                double www = h44->GetBinContent(i, j);
                int ks = 0;

                while (n < 30  && ks <= 5)
                {
                    ks++;
                    if( i > ks)
                    {
                        int ix1 = i -ks;
                        n += h1->GetBinContent(ix1, j);
                        nw += h10->GetBinContent(ix1, j);
                        wl += h3->GetBinContent(ix1 , j);
                        wll += h33->GetBinContent(ix1, j);
                        ww += h4->GetBinContent(ix1, j);
                        www += h44->GetBinContent(ix1, j);
                    }
                    if( i + ks < h1->GetNbinsX())
                    {
                        int ix2 = i+ ks;
                        n += h1->GetBinContent(ix2, j);
                        nw += h10->GetBinContent(ix2, j);
                        wl += h3->GetBinContent(ix2 , j);
                        wll += h33->GetBinContent(ix2, j);
                        ww += h4->GetBinContent(ix2, j);
                        www += h44->GetBinContent(ix2, j);
                    }

                }                
                if( n >= 30)
                {
                    sh1->SetBinContent(i, j, n);
                    sh10->SetBinContent(i, j, nw);
                    sh3->SetBinContent(i, j, wl);
                    sh33->SetBinContent(i, j, wll);
                    sh4->SetBinContent(i, j, ww);
                    sh44->SetBinContent(i, j, www);
                }
                for(int iy = j-3; iy <= j+3; iy++)
                {
                    double w = (iy < j)? 1.-(j - iy)/4. :(iy > j) ? 1-(iy -j)/4. :0.;
                    if(iy >=0 && iy < h1->GetNbinsY())
                    {
                        ns += h1->GetBinContent(i, iy);
                        nws += h10->GetBinContent(i, iy);
                        we += h2->GetBinContent(i, iy);
                        wee += h22->GetBinContent(i,iy);
                    }
                }
                if( ns >= 30)
                {
                    sh11->SetBinContent(i, j, ns);
                    sh101->SetBinContent(i, j, nws);
                    sh2->SetBinContent(i, j, we);
                    sh22->SetBinContent(i, j, wee);
                } 
            }
        }
        for(int i = 0; i < h1->GetXaxis()->GetNbins(); i++)
        {
            for(int j = 0; j < h1->GetYaxis()->GetNbins(); j++)
            {
		    if(sh11->GetBinContent(i, j) > 2 && sh101->GetBinContent(i, j) > 0)
		    {
                    double nw = sh101->GetBinContent(i,j);
                    double ew = sh2->GetBinContent(i,j);
                    he->Fill(sh11->GetXaxis()->GetBinCenter(i), sh11->GetYaxis()->GetBinCenter(j), ew/nw);
                    double eew = sh22->GetBinContent(i,j);
                    hee->Fill(sh11->GetXaxis()->GetBinCenter(i), sh11->GetYaxis()->GetBinCenter(j), sqrt(eew/nw - pow(ew/nw, 2)));
	    	}
            if( sh1->GetBinContent(i, j) > 2 && sh10->GetBinContent(i, j) >0)
            {
                    double n = sh10->GetBinContent(i, j);
                    double lw = sh3->GetBinContent(i, j);
                    double lww = sh33->GetBinContent(i, j);
                    double ww = sh4->GetBinContent(i,j);
                    double www = sh44->GetBinContent(i,j);

                    hl->Fill(sh1->GetXaxis()->GetBinCenter(i), sh1->GetYaxis()->GetBinCenter(j) , lw /n);
                    hll->Fill(sh1->GetXaxis()->GetBinCenter(i), sh1->GetYaxis()->GetBinCenter(j) ,sqrt(lww/n - pow( lw /n, 2)));
                    hw->Fill(sh1->GetXaxis()->GetBinCenter(i), sh1->GetYaxis()->GetBinCenter(j) , ww /n);
                    hww->Fill(sh1->GetXaxis()->GetBinCenter(i), sh1->GetYaxis()->GetBinCenter(j) ,sqrt(www/n - pow( ww /n, 2)));
            }
		
            }
        }
        input_root->Close();

        out_root->cd();
        h1->Write();
        h10->Write();
        h2->Write();
        h22->Write();
        h3->Write();
        h33->Write();
        h4->Write();
        h44->Write();
        sh1->Write();
        sh10->Write();
        sh2->Write();
        sh22->Write();
        sh3->Write();
        sh33->Write();
        sh4->Write();
        sh44->Write();
        he->Write();
        hee->Write();
        hl->Write();
        hll->Write();
        hw->Write();
        hww->Write();
        out_root->Close();



    }

    


    

}
