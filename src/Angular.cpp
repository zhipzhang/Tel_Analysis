#include "TRecData.h"
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TGraph.h"
#include <vector>



int main(int argc, char** argv)
{
    std::string  out_file = "dst.root";
    std::string  in_file = argv[1];
    TRecData* rec = new TRecData();
    TFile* input_root = TFile::Open(in_file.c_str(), "read");
    TTree* rec_tree = (TTree*)input_root->Get("rec_data");

    rec->InitRead(rec_tree);
    TFile* out_root = TFile::Open(out_file.c_str(), "recreate");
    TH2D* h1 = new TH2D("h1", "Direction Error versus Energy",25, -2, 3, 400, 0, 2);
    for(int i = 0; i < rec_tree->GetEntries(); i++)
    {
        h1->Fill(log10(rec->energy), rec->direction_error);
    }
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> sum;
    for(int i = 0; i < h1->GetNbinsX() ;i++)
    {
        double tmp_sum;
        for(int j = 0; j < h1->GetNbinsY(); j++)
        {
            tmp_sum += h1->GetBinContent(i,j);
        }
        sum.push_back(tmp_sum);
        x.push_back(h1->GetXaxis()->GetBinCenter(i));
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
                break;
                tmp3 = j;
            }

        }
        tmp2 = h1->GetYaxis()->GetBinUpEdge(tmp3) - (tmp - sum[i] * 0.68)/h1->GetBinContent(i, tmp3) * h1->GetYaxis()->GetBinWidth(1);
        y.push_back(tmp2);
    }
    TGraph* g1 = new TGraph(x.size(), &x[0], &y[0]);
    h1->Write();
    g1->Write();




}
