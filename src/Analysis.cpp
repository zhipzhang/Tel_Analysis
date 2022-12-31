/*

        Author: Zhipeng Zhang 
        Email : zhipzhang@mail.ustc.edu.cn

        Shower Reconstrution : if not set the lookup file, Only the direction rec and core_position will be made
        Like the read_hess program, you can use the --auto-lookup to automatic use the MakeLookup program to generate the Lookup tables
        ! Remember all data stored are Deg and MM , but what functions calls should be Rad and M!!!
*/




























#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <map>



#include "Limits_defined.h"
#include "LACTree.h"
#include "Detect_config.h"
#include "TGraph.h"
#include "TImage_Parameter.h"
#include "TMath.h"
#include "TMcData.h"
#include "TMultiGraph.h"
#include "TRecData.h"
#include "TUserCuts.h"
#include "moments.h"
#include "TH2Poly.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TPaveText.h"
#include "TH2D.h"
#include "THists.h"
#include "TRandom2.h"
#include "TLine.h"
#include "rec_tools.h"
#ifdef IHEP
std::string prefix = "root://eos01.ihep.ac.cn/";
#else
std::string prefix = "";
#endif

void display(TImage_Parameter* img, double pe[][LACT_MAXPIXELS], std::vector<int> *pixel_in_image, float xpix[][LACT_MAXPIXELS], float ypix[][LACT_MAXPIXELS],double *pix_size, TRecData* rec , int ievent);
using namespace std;
//Image Analysis and Shower Reconstruction

static int pix_shape;
static int image_num = 0;
void syntax();
void compute_pixel_neighbors(map<int, vector<int> > *pixel_neighbor, int ntel, int* npix, float x_pix[][LACT_MAXPIXELS], float y_pix[][LACT_MAXPIXELS],
                          double *pix_size);
void EnergyRec( TRecData* rec, TImage_Parameter* image, TH2D*, TH2D*);
bool Compute_Shape(TRecData* rec, TImage_Parameter* image, TH2D*, TH2D*, TH2D*, TH2D*);
void FillHistogram(TImage_Parameter* image, TRecData* rec);
void InitHistogram();

int main(int argc, char** argv)
{
    std::string input_file = "";
    std::string out_file = "dst.root";
    TImage_Parameter* image = new TImage_Parameter();
    TRandom2* tr3 = new TRandom2();
    TUserCuts *cuts_option = new TUserCuts();
    TRecData* rec = new TRecData();

    TH2D  *he, *hee, *hl, *hll, *hw , *hww;
     

    bool do_image_clean = true; // true if we do the tail-cuts image clean
    bool clean2 = false;
    double auto_lookup = false;
    double erange[3] {1, 10, 100};
    double num_weight[4]{1.0, 1.0, 1.0, 1.0};
    double tail_cuts[2]={5,10} ;

    int ntel;
    static double tel_pos[LACT_MAXTEL][3];
    double focal_length[LACT_MAXTEL];
    double effective_focal_length[LACT_MAXTEL];
    double pe_list[LACT_MAXTEL][LACT_MAXPIXELS]{0};
    int npix[LACT_MAXTEL];
    double tel_direction[LACT_MAXTEL][2];
    static float x_pix[LACT_MAXTEL][LACT_MAXPIXELS];
    static float y_pix[LACT_MAXTEL][LACT_MAXPIXELS];
    double pix_size;                                    // we suppose that all pixel have same pix_size
      //Use Static otherwise the Stack size will exceed

    vector<int> pixel_in_image[LACT_MAXTEL];               //pixel id  which pass the tail-cuts 
    map<int, vector<int> > pixel_neighbors[LACT_MAXTEL];  //use map to store the pixel_neighbors
    bool exist_lookup = false;                           //set it to True will make Tree have the MRSW and MRSL branch
    bool do_shape_cut = false;
    int n_img;
    int n_img2;
    string lookup_file = "";
    string hist_file = "hist.root";

    while(argc > 2)
    {
        if(strcmp(argv[1], "--help") == 0)
        {
            syntax();
            argc--;
            argv++;
            continue;
        }
        if(strcmp(argv[1], "--tail-cuts") == 0)
        {
            int na = 0;
            do_image_clean = 1;
            if( (na = sscanf(argv[2], "%lf,%lf", &tail_cuts[0], &tail_cuts[1]) < 2))
            {
                std::cerr << "Syntax error in tail-cuts parameter" << std::endl;
            }
            argc -= 2;
            argv += 2;
            continue;
        }
        if(strcmp(argv[1], "--lookup-file") == 0)
        {
            exist_lookup = 1;
            lookup_file = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        if(strcmp(argv[1], "--out_file") == 0)
        {
            out_file = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        if(strcmp(argv[1], "--auto-lookup") == 0)
        {
            argc -= 1;
            argv += 1;
            auto_lookup = true;
            continue;
        }
        if(strcmp(argv[1], "--shape-cut") == 0)
        {
            float w_cut[2];
            float l_cuts[2];
            size_t n = sscanf(argv[2], "%f,%f,%f,%f", &w_cut[0], &w_cut[1], &l_cuts[0], &l_cuts[1]);
            if( n < 4)
            {
                fprintf(stderr, "Syntax Error ! you have to specify all four parameters for --shape-cut \n ");
                exit(1);
            }
            do_shape_cut = true;
            cuts_option->SetWidthCut(w_cut);
            cuts_option->SetLengthCut(l_cuts);
            argc -= 2;
            argv += 2;
            continue;
        }
        if(strcmp(argv[1], "--num-weight") == 0)
        {
            size_t n = sscanf(argv[2], "%lf,%lf,%lf,%lf", &num_weight[0], &num_weight[1], &num_weight[2], &num_weight[3]);
            if( n < 4)
            {
                fprintf(stderr, "Syntax Error ! you have to specify all four parameters for --num-weight");
                exit(1);
            }
            argc -= 2;
            argv += 2;
            continue;
        }
        if(strcmp(argv[1], "--histogram-file") == 0)
        {
            hist_file = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        if(strcmp(argv[1], "--image-clean") == 0)
        {
            clean2 = true;
            do_image_clean = false;
            std::cout << "clean2 will override  the tail-cuts clean" << std::endl;
            argc -= 1;
            argv += 1;
            continue;
        }
        else
        {
            break;
        }
    }
    if(clean2)
    {
        do_image_clean = false;
    }
    //  Try to Open the lookupfile
    if(!lookup_file.empty() && exist_lookup)
    {
        std::cout << "Tring to Open lookupfile " << lookup_file << std::endl;
        TFile* lookup_root = TFile::Open(lookup_file.c_str(), "read");
        if(lookup_root->IsZombie())
        {
            syntax();
            perror(lookup_file.c_str());
            exit(EXIT_FAILURE);
        }
    rec->SetLookup();
    he = (TH2D*) ( lookup_root->Get("he"));
    hee = (TH2D*)lookup_root->Get("hee");
    hl = (TH2D*)lookup_root->Get("hl");
    hw = (TH2D*)lookup_root->Get("hw");
    hll = (TH2D*)lookup_root->Get("hll");
    hww = (TH2D*)lookup_root->Get("hww");
    }

    TFile *out_root = TFile::Open(out_file.c_str(), "recreate");
    //TH1D* h4 = new TH1D("h4", "True Energy Histogram (After Clean and Pass Size Cut events) In 400m Radius", 20, -1, 3);
    //TH1D* h5 = new TH1D("h5", "True Energy Histogram (After Clean and Pass Size Cut events) In 450m Radius", 20, -1, 3);
    TTree *image_tree = new TTree("image_parameter", "image_parameter");
    image_tree->Branch("image_parameter",&image );
    rec->InitWrite();
    
    THists *hists = new THists(hist_file);

    while( argc > 1 || !input_file.empty() )
    {
        std::cout<< "Opening the input file" <<argv[1] << std::endl;;
        input_file = prefix + argv[1];
        argc--;
        argv++;
        TFile* root_file = TFile::Open(input_file.c_str(), "read");
        TTree* event_tree = (TTree*) root_file->Get("event");
        TTree* mc_tree = (TTree*) root_file->Get("mc");
        TTree* config_tree = (TTree*) root_file->Get("config_tree");
        input_file = "";

        float energy, weight;
        mc_tree->SetBranchAddress("MCe0", &energy);
        for(int i = 0; i < mc_tree->GetEntries(); i++)
        {
            mc_tree->GetEntry(i);
            if(energy < 1 )
                weight = pow(energy, -1.7) * num_weight[0];
            else if (energy < 10)
                weight = pow(energy, -1.7) * num_weight[1];
            else if (energy < 100)
                weight = pow(energy, -1.7) * num_weight[2];
            else
                weight = pow(energy, -1.7) * num_weight[3];
            hists->h1->Fill(log10(energy), weight);

        }
        LACTree* DSTTree = new LACTree();
        Detect* detect = new Detect();
        
        DSTTree->initEventTree(event_tree);
        detect->InitRead(config_tree);

        if(config_tree )        
        {
            for(int i = 0; i < config_tree->GetEntries(); i++)
            {
                config_tree->GetEntry(i);
                if(detect->GetNtel() > LACT_MAXTEL)
                {
                    cerr << "There are too many tels When configuring" << endl;
                    exit(EXIT_FAILURE);
                }
                if(detect->GetNpix() > LACT_MAXPIXELS)
                {
                    cerr << "There are too many pixels when configuring " << endl;
                    exit(EXIT_FAILURE);
                }
                ntel = detect->GetNtel();
                int itel = detect->GetTel_id() - 1;
                tel_pos[i][0] = detect->Pos[0];
                tel_pos[i][1] = detect->Pos[1];
                tel_pos[i][2] = detect->Pos[2];
                focal_length[itel ] = detect->GetFocal_Length();
                effective_focal_length[itel ] = detect->GetEffective_Focal_length();
                npix[itel ] = detect->GetNpix();
                for(int p = 0; p < detect->GetNpix(); p++)
                {
                    x_pix[itel][p] = detect->X_PixMM[p] / 1000 ;
                    y_pix[itel][p] = detect->Y_PixMM[p] / 1000 ;
                }
                pix_size = detect->Pixel_Size[0] /1000;
                pix_shape = detect->Pixel_Shape[0];


            }
           compute_pixel_neighbors(pixel_neighbors, ntel, npix, x_pix, y_pix, &pix_size);

            
        }
        else
        {
            cerr << "No config Tree in the ROOT File" << input_file << endl;
        }

        // Get Telescope position 
        rec->SetTelPos(tel_pos);
        int ntrig;
        for( int i = 0 ; i < event_tree->GetEntries(); i++)
        {
            bool shape_ok = false;
            event_tree->GetEntry(i);
            rec->GetData(DSTTree);
            rec->ReWeight(erange, num_weight);
            ntrig = DSTTree->GetNtrig();
            if(!DSTTree->fadc_read_write )
            {
                if( i == 0)
                    std::cout << "Be careful ! We don't have ADC data so we will use ideal Pe data instead" << std::endl;
                if(!DSTTree->fillPeLeaf)
                {
                    std::cerr << "Can not Find Pe data !" << std::endl;
                    //exit(EXIT_FAILURE);
                }
                else
                {
                   for(int j = 0; j < DSTTree->GetNtrig(); j++)
                   {
                        int tel_id = DSTTree->GetTrigI(j) - 1;
                        for(int p = 0; p < npix[tel_id]; p++)
                        {
                            pe_list[tel_id][p] =  DSTTree->pe_list[j][p] + tr3->Poisson(8.91) - 8.91   ;
                            if(pe_list[tel_id][p] < 0 )
                            {
                                pe_list[tel_id][p] = 0;
                            }
                        }
                        tel_direction[tel_id][0] = DSTTree->Point_Az[j];
                        tel_direction[tel_id][1] = DSTTree->Point_Al[j];
                        
                   } 
                    
                }
            }
            else //if we have the fadc data, we can get pe from it
            {
                
            }
            rec->SetTelDirection(tel_direction);
            hists->h2->Fill(log10(rec->GetEnergy()), rec->weight);
            if(do_image_clean)
            {
                image_clean(pe_list, pixel_neighbors, tail_cuts, pixel_in_image, DSTTree->GetNtrig(), DSTTree->GetTrig_List(), npix);
            }
            
            else {
                image_clean2(pe_list, pixel_neighbors, tail_cuts, pixel_in_image, DSTTree->GetNtrig(), DSTTree->GetTrig_List(), npix);
            }
            compute_moments(image, rec, pe_list, pixel_in_image, ntrig, DSTTree->GetTrig_List(), x_pix, y_pix);
            image->ConvertToRad(focal_length);
             
             
            
            if(cuts_option->Check(image))
            {
                hists->h3->Fill(log10(rec->GetEnergy()), rec->weight);
                hists->h4->Fill(rec->GetCoredist(), log10(rec->GetEnergy()), rec->weight);
            }
            else
            {
                memset(pe_list, 0, LACT_MAXTEL*LACT_MAXPIXELS*sizeof(pe_list[0][0]));
                for(int i = 0; i< LACT_MAXTEL; i++)
                {
                    pixel_in_image[i].clear();
                }
                rec->Reset();
                image->clear();
                continue;
            }
            if(!rec->RecShower(image, cuts_option))
            {
                memset(pe_list, 0, LACT_MAXTEL*LACT_MAXPIXELS*sizeof(pe_list[0][0]));
                for(int i = 0; i< LACT_MAXTEL; i++)
                {
                    pixel_in_image[i].clear();
                }
                rec->Reset();
                image->clear();
                continue;
            }
            ComputeRp(image, rec, tel_pos);
           
            
           
            
            rec->compute_direction_error();
            /*
            for(int h = 0 ; h < image->image_tel.size(); h++)
            {
                int tel_id = image->image_tel[h];
                if( image->GetDist(tel_id) * TMath::RadToDeg() >4 && image_num < 100)
                {
                    display(image, pe_list, pixel_in_image, x_pix, y_pix, &pix_size, rec, i);
                    image_num++;
                    break;

                }
            }
            */
            /*
            if( i == 2506)
            {
                std::cout << "image x is "<< image->GetTelImageX(1) << std::endl;
                std::cout << "image y is "<< image->GetTelImageY(1) << std::endl;
                std::cout << "alpha is "<< image->GetTelAlpha(1) << std::endl;
                std::cout << "length is "<< image->GetTelLength(1) << std::endl;
                std::cout << "width is "<< image->GetTelwidth(1) << std::endl;
                std::cout << rec->GetDirectionError()<< "direction_error";
                display(image, pe_list, pixel_in_image, x_pix, y_pix, &pix_size, rec, i);

            }
            */
            hists->h60->Fill(rec->GetDirectionError(), rec->weight);
            hists->h62->Fill(log10(rec->energy), rec->direction_error, rec->weight);
            int ie = (int)(10. * log10(rec->GetEnergy()/0.01)/3.) - 4 ;
            if(ie > 10)
            {
                ie = 10;
            }

            if(ie >= 0)
            {
                hists->h100[ie]->Fill(rec->GetCorex(), rec->GetCorey(), rec->weight);
                for(int i = 0; i < image->image_tel.size(); i++)
                {
                    int tel_id = image->image_tel[i];
                    double phi = atan((image->GetTelImageY(tel_id) - image->source_x[tel_id] ) / (image->GetTelImageY(tel_id) - image->source_y[tel_id]));
                    double miss = line_point_distance(image->GetTelImageX(tel_id), image->GetTelImageY(tel_id), 0, cos(image->GetTelAlpha(tel_id)), sin(image->GetTelAlpha(tel_id)), 0, image->source_x[tel_id], image->source_y[tel_id], 0);
                    hists->h200[ie]->Fill(image->GetTelRp(tel_id), image->GetDist(tel_id) * TMath::RadToDeg(), rec->weight);
                    hists->h250[ie]->Fill(image->GetDist(tel_id) * TMath::RadToDeg(), fabs(phi - image->GetTelAlpha(tel_id)) * TMath::RadToDeg());
                    hists->h2250[ie]->Fill(image->GetDist(tel_id) * TMath::RadToDeg(), miss * TMath::RadToDeg());
                    hists->h450[ie]->Fill(image->GetTelRp(tel_id), image->GetTelwidth(tel_id) * TMath::RadToDeg());
                    hists->h550[ie]->Fill(image->GetTelRp(tel_id), image->GetTelLength(tel_id) * TMath::RadToDeg());
                    hists->hist_max->Fill(rec->xmax, image->GetTelRp(tel_id), image->GetDist(tel_id) * TMath::RadToDeg());
                }
                hists->h600[ie]->Fill(rec->GetCoredist(), rec->GetDirectionError(), rec->weight);

            } 
            if(exist_lookup)
            {
                EnergyRec(rec, image, he, hee);
                bool rc = Compute_Shape(rec, image, hw, hww, hl, hll);
                if(!rc)
                {
                    memset(pe_list, 0, LACT_MAXTEL*LACT_MAXPIXELS*sizeof(pe_list[0][0]));
                    for(int i = 0; i< LACT_MAXTEL; i++)
                    {
                        pixel_in_image[i].clear();
                    }
                    rec->Reset();
                    image->clear();
                    
                    continue;

                }
                if(do_shape_cut &&rec->MRSL >cuts_option->length_cut[0] && rec->MRSL < cuts_option->length_cut[1] 
                   && rec->MRSW >cuts_option->width_cut[0] && rec->MRSW <cuts_option->width_cut[1])
                   {
                    shape_ok = true;
                   }
                hists->h303->Fill(rec->MRSW, rec->weight);
                hists->h304->Fill(rec->MRSL, rec->weight);
                hists->h301->Fill(rec->MRSW, rec->MRSL, rec->weight);
                if(rec->direction_error  < 1.0)
                {
                    hists->h302->Fill(rec->MRSW, rec->MRSL, rec->weight);
                }
            }
            if(ie >= 0)
            {
                if(exist_lookup)
                {
                    hists->h400[ie]->Fill(rec->MRSW, rec->weight);
                    hists->h500[ie]->Fill(rec->MRSL, rec->weight);
                }
            }
            if(shape_ok)
            {
                hists->h5->Fill(log10(rec->energy), rec->weight);
                hists->h63->Fill(log10(rec->energy), rec->direction_error, rec->weight);
                hists->h61->Fill(rec->direction_error, rec->weight);
                hists->h64->Fill(rec->theta2, rec->weight);
                if(rec->theta2 < 0.0125)
                {
                    hists->h6->Fill(log10(rec->energy), rec->weight);
                }
            }
            image_tree->Fill();
            rec->GetRecTree()->Fill();
            memset(pe_list, 0, LACT_MAXTEL*LACT_MAXPIXELS*sizeof(pe_list[0][0]));
            for(int i = 0; i< LACT_MAXTEL; i++)
            {
                pixel_in_image[i].clear();
            }
            rec->Reset();
            image->clear();
        }
            delete DSTTree;
            delete detect;
            root_file->Close();
            for(int i = 0 ; i < LACT_MAXTEL; i++)
            {
                pixel_neighbors[i].clear();
            }

    }
    hists->Write();
    out_root->cd();
    image_tree->Write();
    rec->GetRecTree()->Write();
    out_root->Write();
    out_root->Close();

}
void syntax()
{
    printf("--help show help message \n");
    printf("--tail-cuts  use two level cuts to clean image \n");

}
void compute_pixel_neighbors(map<int, vector<int> > *pixel_neighbor, int ntel, int* npix, float x_pix[][LACT_MAXPIXELS], float y_pix[][LACT_MAXPIXELS], 
                            double* pix_size)
{
    double x_j, y_j;
    double x_k, y_k;
    double size_j, size_k;
    for(int itel = 0; itel < ntel; itel++)
    {
        for(int j = 0; j < npix[itel]; j++)
        {
            x_j = x_pix[itel][j];
            y_j = y_pix[itel][j];
            size_j = *pix_size * 0.5;
            for(int k = 0; k < npix[itel]; k++)
            {
                if( k != j)
                {
                    x_k = x_pix[itel][k];
                    y_k = y_pix[itel][k];
                    size_k = *pix_size * 0.5;

                    if(sqrt(pow(x_j - x_k, 2) + pow(y_j - y_k ,2)) < 1.5 *(size_j +size_k) )
                    {
                        if(pixel_neighbor[itel][j].size() < MAX_NEIGHBOR)
                            pixel_neighbor[itel][j].push_back(k);
                        else
                        {
                            break;
                        }
                    }

                }
            }
            
        }
    }
}
void display(TImage_Parameter* image, double pe[][LACT_MAXPIXELS], std::vector<int> *pixel_in_image, float xpix[][LACT_MAXPIXELS], float ypix[][LACT_MAXPIXELS], double *pix_size, TRecData* rec, int ievent)
{
    int num = 0;
    for(int i = 0; i < image->image_tel.size(); i++)
    {
        int tel_id = image->image_tel[i];
        if(image->GetTelSize(tel_id) > 50)
        {
            double source_x, source_y;
            angles_to_offset(rec->GetAzimuth() * TMath::DegToRad(), rec->GetAltitude() * TMath::DegToRad(), rec->Tel_direction[tel_id][0] * TMath::DegToRad(), rec->Tel_direction[tel_id][1] * TMath::DegToRad()
                , image->GetTelFocal(tel_id), &source_x, &source_y);
            source_x = source_x/image->GetTelFocal(tel_id) * TMath::RadToDeg();
            source_y = source_y/image->GetTelFocal(tel_id) * TMath::RadToDeg();

            
            
            TCanvas* camera_image = new TCanvas(Form("LACT image of camera %d", tel_id), "LACT image", 1600, 1600);
            TH2Poly* camera = new TH2Poly(Form("camera %d", tel_id), "", -6, 6, -6, 6);
            if(pix_shape != 1)
            {
                for(auto j: pixel_in_image[tel_id])
                {
                    double binsize = *pix_size / image->GetTelFocal(tel_id) * TMath::RadToDeg();
                    double x = xpix[tel_id][j]/ image->GetTelFocal(tel_id) * TMath::RadToDeg();
                    double y = ypix[tel_id][j]/ image->GetTelFocal(tel_id) * TMath::RadToDeg();
                    double bin_x[4] = {x - 0.5 * binsize, x + 0.5*binsize, x + 0.5*binsize, x - 0.5*binsize};
                    double bin_y[4] = {y - 0.5 * binsize, y - 0.5*binsize, y + 0.5*binsize, y + 0.5*binsize};
                    camera->AddBin(4, bin_x, bin_y);
                    camera->Fill(x, y, pe[tel_id][j]);
                }
            }
            else
            {
                for( auto j : pixel_in_image[tel_id])
                {
                    double binsize = *pix_size / image->GetTelFocal(tel_id) * TMath::RadToDeg() * 0.5 ;
                    double x = xpix[tel_id][j]/ image->GetTelFocal(tel_id) * TMath::RadToDeg();
                    double y = ypix[tel_id][j]/ image->GetTelFocal(tel_id) * TMath::RadToDeg();
                    double bin_x[6] = {x - binsize, x - 0.5 * binsize , x + 0.5*binsize , x + binsize, x+0.5*binsize, x -0.5*binsize};
                    double bin_y[6] = {y, y - 0.5*binsize*TMath::Sqrt(3), y - 0.5*binsize*TMath::Sqrt(3) , y, y+0.5*binsize*TMath::Sqrt(3), y+ 0.5*binsize*TMath::Sqrt(3)};
                    camera->AddBin(6, bin_x, bin_y);
                    camera->Fill(x, y, pe[tel_id][j]);
                }

            }
            camera->SetStats(0);
             
            camera->Draw();

            TGraph* g1 = new TGraph();
            g1->AddPoint(source_x, source_y);
            g1->SetMarkerStyle(4);
            g1->SetMarkerSize(4);
            g1->SetMarkerColor(2);
            TGraph* g2 = new TGraph();
            g2->AddPoint(rec->rec_x/image->GetTelFocal(tel_id) * TMath::RadToDeg(), rec->rec_y/image->GetTelFocal(tel_id) * TMath::RadToDeg());
            g2->SetMarkerStyle(5);
            g2->SetMarkerColor(2);
            g2->SetMarkerSize(4);

            TMultiGraph* mg = new TMultiGraph();
            mg->Add(g1);
            mg->Add(g2);
            camera_image->cd();
            mg->Draw("p");

            TEllipse* ellipse = new TEllipse(image->GetTelImageX(tel_id) * TMath::RadToDeg(), image->GetTelImageY(tel_id) * TMath::RadToDeg(), image->GetTelLength(tel_id) * TMath::RadToDeg(), image->GetTelwidth(tel_id) * TMath::RadToDeg(),
                                            0, 360, image->GetTelAlpha(tel_id)  * TMath::RadToDeg() );
            ellipse->SetLineWidth(2);
            ellipse->SetLineColor(2);
            ellipse->SetFillStyle(0);
            ellipse->Draw();
            
            TLine* line = new TLine(image->GetTelImageX(tel_id) * TMath::RadToDeg() - 5 * TMath::Cos(image->GetTelAlpha(tel_id)), image->GetTelImageY(tel_id) * TMath::RadToDeg() - 5 * TMath::Sin(image->GetTelAlpha(tel_id)), image->GetTelImageX(tel_id) * TMath::RadToDeg() + 5 * TMath::Cos(image->GetTelAlpha(tel_id)), image->GetTelImageY(tel_id) * TMath::RadToDeg() + 5 * TMath::Sin(image->GetTelAlpha(tel_id)));
            line->SetLineWidth(2);
            line->SetLineColor(4);
            line->Draw();

            
            TEllipse* el2 = new TEllipse(0, 0, 5.0, 5.0, 0, 360);
            el2->SetLineWidth(2);
            el2->SetLineColor(kBlack);
            el2->SetFillStyle(0);
            el2->Draw();
            
            TPaveText *pavet = new TPaveText(-6, 6.3, 6, 7.6);
            pavet->SetFillStyle(0);
            pavet->AddText(Form("event_number: %d, Tel: %d, energy: %.4lf , azimuth: %.4lf, zenith: %.4lf Rp:%.4lf m size: %f core_pos : %fm ,%fm  dist %f hmax %lf emax %lf", 
                  rec->GetEventNumber(),tel_id, rec->GetEnergy(), rec->GetAzimuth(),  90 - rec->GetAltitude(), image->GetTelRp(tel_id), image->GetTelSize(tel_id), rec->core_pos[0], rec->core_pos[1], image->GetDist(tel_id) * TMath::RadToDeg() , rec->hmax, rec->emax));
            pavet->Draw("same");
            camera_image->SaveAs(Form("./images4/%dimage_camera%d.png", ievent,tel_id));
        }

    }
}

void EnergyRec(TRecData* rec, TImage_Parameter* image, TH2D* h1, TH2D* h2)
{
    double all_e = 0.;
    double all_w = 0.;
    for(int i = 0 ; i< image->nsig; i++)
    {
        int itel = image->sig_tel[i];
        double size = image->GetTelSize(itel);
        double rp = image->GetRecRp(itel);
        int xbin = h1->GetXaxis()->FindBin(rp);
        int ybin = h1->GetYaxis()->FindBin(log10(size));
        if(xbin <=0 || xbin > h1->GetNbinsX() || ybin <= 0 || ybin > h1->GetNbinsY())
        {
            continue;
        }
        double energy = size/h1->GetBinContent(xbin, ybin) ;
        double var_e;
        if(h1->GetBinContent(xbin, ybin) > 0)
            var_e =  h2->GetBinContent(xbin, ybin)/h1->GetBinContent(xbin, ybin);
        else
            var_e = 999.;
        double weight = 1/(0.01 + var_e * var_e);
        all_e += log(energy) * weight;
        all_w += weight;
    }
    if(all_w == 0)
        rec->SetRecEnergy(-1);
    else
        rec->SetRecEnergy(exp(all_e/all_w));
}

//function to Compute MRSW MRSL
bool Compute_Shape(TRecData* rec, TImage_Parameter* image, TH2D* hw, TH2D* hww, TH2D* hl, TH2D* hll)
{
    double all_rsw = 0.,all_rsl = 0.;
    int n = image->nsig;
    for( int i = 0; i < image->nsig  ; i++)
    {
        int tel_id = image->sig_tel[i];
        double size = image->GetTelSize(tel_id);
        double rp = image->GetRecRp(tel_id);
        int xbin = hw->GetXaxis()->FindBin(rp);
        int ybin = hw->GetYaxis()->FindBin(log10(size));
        if( xbin <= 0 || xbin >=hw->GetNbinsX() || ybin <=0 || ybin > hw->GetNbinsY() || hww->GetBinContent(xbin,ybin) < 1e-8 || hll->GetBinContent(xbin, ybin) < 1e-8)
        {
            n--;
            continue;
        }
        all_rsw += (image->GetTelwidth(tel_id) - hw->GetBinContent(xbin, ybin)) / hww->GetBinContent(xbin, ybin);
        all_rsl += (image->GetTelLength(tel_id) - hl->GetBinContent(xbin, ybin)) / hll->GetBinContent(xbin, ybin);
    }
    if( n > 0)
    {
        rec->SetShape(all_rsw/n , all_rsl/n);
        return true;
    }
    else
    {
        return false;
    }


}
