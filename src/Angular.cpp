#include "TRecData.h"
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TGraph.h"
#include <vector>
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"



int main(int argc, char** argv)
{
    std::string  out_file = "dst.root";
    std::string  in_file = argv[1];
    TRecData* rec = new TRecData();
    TFile* input_root = TFile::Open(in_file.c_str(), "read");
    TTree* rec_tree = (TTree*)input_root->Get("rec_data");

    rec->InitRead(rec_tree);
    TFile* out_root = TFile::Open(out_file.c_str(), "recreate");
    out_root->cd();
    TH2D* h1 = new TH2D("h1", "Direction Error versus Energy", 20, -0.8, 3, 4000, 0, 10);
    TH2D* h2 = new TH2D("h2", "Direction Error versus Energy", 20, -0.8, 3, 4000, 0, 10);
    for(int i = 0; i < rec_tree->GetEntries(); i++)
    {
	rec_tree->GetEntry(i);
        if(rec->GetCoredist() < 250)
		h2->Fill(log10(rec->energy), rec->direction_error);
        h1->Fill(log10(rec->energy), rec->direction_error);
    }
    std::vector<double> x;
    std::vector<double> y; 
    std::vector<double> sum; 
    for(int i = 0; i < h1->GetNbinsX(); i++)
    {
	double tmp_sum = 0;
        for(int j = 0; j < h1->GetNbinsY(); j++)
            tmp_sum +=h1->GetBinContent(i, j);
        sum.push_back(tmp_sum);
        x.push_back(h1->GetXaxis()->GetBinCenter(i));
        std::cout << x[i]<< std::endl;
    }
    for(int i = 0; i < h1->GetXaxis()->GetNbins(); i++)
    {
        
        double tmp = 0;
        double tmp2 = 0;
        int tmp3 = 0;
        for(int j = 0; j < h1->GetYaxis()->GetNbins(); j++)
        {
            tmp += h1->GetBinContent(i,j);
            if( tmp > sum[i] * 0.68)
            {
                tmp3 = j;
                break;
            }

        }
        tmp2 = h1->GetYaxis()->GetBinUpEdge(tmp3) - (tmp - sum[i] * 0.68)/h1->GetBinContent(i, tmp3) * h1->GetYaxis()->GetBinWidth(1);
        y.push_back(tmp2);
        std::cout << y[i]<< std::endl;
    }
    std::vector<double> x2;
    std::vector<double> y2; 
    std::vector<double> sum2; 
    for(int i = 0; i < h2->GetNbinsX(); i++)
    {
	double tmp_sum = 0;
        for(int j = 0; j < h2->GetNbinsY(); j++)
            tmp_sum +=h2->GetBinContent(i, j);
        sum2.push_back(tmp_sum);
        x2.push_back(h1->GetXaxis()->GetBinCenter(i));
        std::cout << x2[i]<< std::endl;
    }
    for(int i = 0; i < h2->GetXaxis()->GetNbins(); i++)
    {
        
        double tmp = 0;
        double tmp2 = 0;
        int tmp3 = 0;
        for(int j = 0; j < h2->GetYaxis()->GetNbins(); j++)
        {
            tmp += h2->GetBinContent(i,j);
            if( tmp > sum2[i] * 0.68)
            {
                tmp3 = j;
                break;
            }

        }
        tmp2 = h2->GetYaxis()->GetBinUpEdge(tmp3) - (tmp - sum2[i] * 0.68)/h2->GetBinContent(i, tmp3) * h2->GetYaxis()->GetBinWidth(1);
        y2.push_back(tmp2);
        std::cout << y2[i]<< std::endl;
    }
    double c[15] = {-0.25, -0.05, 0.15, 0.35, 0.55, 0.75, 0.95, 1.15, 1.35, 1.55, 1.75, 1.95, 2.15, 2.35, 2.55};
    double d[15] = {0.086, 0.062, 0.048, 0.0400, 0.035954, 0.032515, 0.029608, 0.028704, 0.027684, 0.027569,0.027110, 0.026549, 0.023432, 0.020169, 0.019838};
    out_root->cd();
    TCanvas *c1 = new TCanvas("c1","angular resolution");
    TGraph* g1 = new TGraph(x.size(), &x[0], &y[0]);
    g1->SetLineColor(kRed);
    TGraph* g3 = new TGraph(x2.size(), &x2[0], &y2[0]);
    auto mg = new TMultiGraph("mg","mg");
    mg->Add(g3);
    mg->Add(g1);
    mg->SetTitle("Angular Resolution;log(Energy/1TeV);Direction Error ");
    mg->Draw("alp");
    TLegend* leg = new TLegend(0.8,0.8,0.95,0.95);
    leg->AddEntry(g1, "All Events");
    leg->AddEntry(g3, "Pass Shape Cuts");
    leg->Draw();
    
    
    TGraph* g2 = new TGraph(15, c, d);
    h1->Write();
    h2->Write();
    c1->Write();
    out_root->Write();
    input_root->Close();
    




}
