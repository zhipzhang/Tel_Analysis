/*

        Author: Zhipeng Zhang 
        Email : zhipzhang@mail.ustc.edu.cn

        Shower Reconstrution : if not set the lookup file, Only the direction rec and core_position will be made
        Like the read_hess program, you can use the --auto-lookup to automatic use the MakeLookup program to generate the Lookup tables
        ! Remember all data stored are Deg and MM , but what functions calls should be Rad and M!!!
*/




























#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>
#include <map>



#include "Limits_defined.h"
#include "LACTree.h"
#include "Detect_config.h"
#include "TImage_Parameter.h"
#include "TMcData.h"
#include "TRecData.h"
#include "TCuts.h"
#include "moments.h"

using namespace std;
//Image Analysis and Shower Reconstruction

void syntax();
void compute_pixel_neighbors(map<int, vector<int> > *pixel_neighbor, int ntel, int* npix, double x_pix[][LACT_MAXPIXELS], double y_pix[][LACT_MAXPIXELS],
                          double (*pix_size)[LACT_MAXPIXELS]);

static LACTree* DSTTree;
static Detect* detect;
int main(int argc, char** argv)
{
    const char* input_file = NULL;
    string out_file = "dst.root";
    TImage_Parameter* image = new TImage_Parameter();
    TCuts *cuts_option = new TCuts();
    TRecData* rec = new TRecData();

    bool do_image_clean = true; // true if we do the tail-cuts image clean
    double auto_lookup = false;
    double tail_cuts[2]{0} ;

    int ntel;
    static double tel_pos[LACT_MAXTEL][3];
    double focal_length[LACT_MAXTEL];
    double effective_focal_length[LACT_MAXTEL];
    static double pe_list[LACT_MAXTEL][LACT_MAXPIXELS]{0};
    int npix[LACT_MAXTEL];
    double tel_direction[LACT_MAXTEL][2];
    static double x_pix[LACT_MAXTEL][LACT_MAXPIXELS];
    static double y_pix[LACT_MAXTEL][LACT_MAXPIXELS];
    static double pix_size[LACT_MAXTEL][LACT_MAXPIXELS];     //Use Static otherwise the Stack size will exceed

    vector<int> pixel_in_image[LACT_MAXTEL];               //pixel id  which pass the tail-cuts 
    map<int, vector<int> > pixel_neighbors[LACT_MAXTEL];  //use map to store the pixel_neighbors
    bool exist_lookup = false;                           //set it to True will make Tree have the MRSW and MRSL branch
    string lookup_file = "";

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
            argc -= (na+1);
            argv += (na+1);
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
        }
        else
        {
            break;
        }
    }

    //  Try to Open the lookupfile
    if(!lookup_file.empty())
    {
        std::cout << "Tring to Open lookupfile " << lookup_file << std::endl;
        TFile* lookup_root = TFile::Open(lookup_file.c_str(), "read");
        if(lookup_root->IsZombie())
        {
            syntax();
            perror(lookup_file.c_str());
            exit(EXIT_FAILURE);
        }
    }

    TFile *out_root = TFile::Open(out_file.c_str(), "recreate");
    TTree *image_tree = new TTree("image_parameter", "image_parameter");
    image_tree->Branch("image_prameter",&image );
    rec->InitWrite();
    

    while( argc > 1 || input_file != NULL)
    {
        input_file = argv[1];
        argc--;
        argv++;
        TFile* root_file = TFile::Open(input_file, "read");
        TTree* event_tree = (TTree*) root_file->Get("event");
        TTree* config_tree = (TTree*) root_file->Get("config_tree");
        input_file = NULL;

        DSTTree = new LACTree();
        detect = new Detect();
        
        DSTTree->initEventTree(event_tree);
        detect->InitRead(config_tree);

        if(config_tree)        
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
                    x_pix[itel][p] = detect->X_PixMM[p] / 1000;
                    y_pix[itel][p] = detect->Y_PixMM[p] / 1000;
                    pix_size[itel][p] = detect->Pixel_Size[p] /1000;
                }


            }
            compute_pixel_neighbors(pixel_neighbors, ntel, npix, x_pix, y_pix, pix_size);
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
            event_tree->GetEntry(i);
            rec->GetData(DSTTree);
            ntrig = DSTTree->GetNtrig();
            if(!DSTTree->fadc_read_write)
            {
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
                            pe_list[tel_id][p] = DSTTree->pe_list[j][p];
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
            if(do_image_clean)
            {
                image_clean(pe_list, pixel_neighbors, tail_cuts, pixel_in_image, DSTTree->GetNtrig(), DSTTree->GetTrig_List(), npix);
            }
            compute_moments(image, pe_list, pixel_in_image, ntrig, DSTTree->GetTrig_List(), x_pix, y_pix);
            image->ConvertToRad(focal_length);
            rec->RecShower(image, cuts_option);
            image_tree->Fill();
            rec->GetRecTree()->Fill();
            printf("%lf %lf \n",rec->rec_core[0], rec->rec_core[1]);
            rec->Reset();
            image->clear();
            
        }

    }
    out_root->cd();
    image_tree->Write();
    rec->GetRecTree()->Write();
    out_root->Close();
}
void syntax()
{
    printf("--help show help message \n");
    printf("--tail-cuts  use two level cuts to clean image \n");

}
void compute_pixel_neighbors(map<int, vector<int> > *pixel_neighbor, int ntel, int* npix, double x_pix[][LACT_MAXPIXELS], double y_pix[][LACT_MAXPIXELS], 
                            double (*pix_size)[LACT_MAXPIXELS])
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
            size_j = pix_size[itel][j] * 0.5;
            for(int k = 0; k < npix[itel]; k++)
            {
                if( k != j)
                {
                    x_k = x_pix[itel][k];
                    y_k = y_pix[itel][k];
                    size_k = pix_size[itel][k] * 0.5;

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