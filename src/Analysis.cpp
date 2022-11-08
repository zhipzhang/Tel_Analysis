#include "TFile.h"
#include "TTree.h"
#include <string>
#include <map>


#include "LACTree.h"
#include "Detect_config.h"
#include "Limits_defined.h"
#include "TImage_Parameter.h"
#include "TRecData.h"
#include "TCuts.h"
#include "moments.h"

using namespace std;
//Image Analysis and Shower Reconstruction
void syntax();
void syntax()
{

}

void compute_pixel_neighbors(map<int, vector<int> > *pixel_neighbor, int ntel, int* npix, double x_pix[][LACT_MAXPIXELS], double y_pix[][LACT_MAXPIXELS],
                            double (*pix_size)[LACT_MAXPIXELS]);
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
            size_j = pix_size[itel][j] * 0.05;
            for(int k = 0; k < npix[itel]; k++)
            {
                if( k != j)
                {
                    x_k = x_pix[itel][k];
                    y_k = y_pix[itel][k];
                    size_k = pix_size[itel][k] * 0.05;

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

int main(int argc, char** argv)
{
    const char* input_file = NULL;
    string out_file = "dst.root";
    TImage_Parameter* image = new TImage_Parameter();
    TCuts *cuts_option = new TCuts();
    bool do_image_clean = false; // true if we do the tail-cuts image clean
    double tail_cuts[2];

    int ntel;
    double tel_pos[LACT_MAXTEL][3];
    double focal_length[LACT_MAXTEL];
    double effective_focal_length[LACT_MAXTEL];
    double pe_list[LACT_MAXTEL][LACT_MAXPIXELS]{0};
    int npix[LACT_MAXTEL];
    double x_pix[LACT_MAXTEL][LACT_MAXPIXELS]{0};
    double y_pix[LACT_MAXTEL][LACT_MAXPIXELS]{0};
    double pix_size[LACT_MAXTEL][LACT_MAXPIXELS]{0};

    vector<int> pixel_in_image[LACT_MAXTEL]; //pixel id  which pass the tail-cuts 
    map<int, vector<int> > pixel_neighbors[LACT_MAXTEL];
    bool exist_lookup = false; //set it to True will make Tree have the MRSW and MRSL branch
    string lookup_file = "";

    while(argc > 2)
    {
        if(strcmp(argv[1], "-h") == 0)
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
        else
        {
            break;
        }
    }
    TFile *out_root = TFile::Open(out_file.c_str(), "recreate");
    TTree *image_tree = new TTree("image_parameter", "image_parameter");
    image_tree->Branch("image_prameter",&image);
    



    while( argc > 1 || input_file != NULL)
    {
        input_file = argv[1];
        argc--;
        argv++;
        TFile* root_file = TFile::Open(input_file, "read");
        TTree* event_tree = (TTree*) root_file->Get("event");
        TTree* config_tree = (TTree*) root_file->Get("config_tree");

        LACTree* DSTTree = new LACTree();
        Detect* detect = new Detect();
        
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
                int itel = detect->GetTel_id();
                focal_length[itel] = detect->GetFocal_Length();
                effective_focal_length[itel] = detect->GetEffective_Focal_length();
                npix[itel] = detect->GetNpix();
                for(int p = 0; p < detect->GetNpix(); p++)
                {
                    x_pix[itel][p] = detect->X_PixMM[p];
                    y_pix[itel][p] = detect->Y_PixMM[p];
                    pix_size[itel][p] = detect->Pixel_Size[p];
                }


            }
            compute_pixel_neighbors(pixel_neighbors, ntel, npix, x_pix, y_pix, pix_size);
        }
        else
        {
            cerr << "No config Tree in the ROOT File" << input_file << endl;
        }
        for( int i = 0 ; i < event_tree->GetEntries(); i++)
        {
            event_tree->GetEntry(i);
            TRecData *rec = new TRecData(DSTTree);
            rec->InitWrite();
            if(!DSTTree->fadc_read_write)
            {
                std::cout << "Be careful ! We don't have ADC data so we will use ideal Pe data instead" << std::endl;
                if(!DSTTree->fillPeLeaf)
                {
                    std::cerr << "Can not Find Pe data !" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else
                {
                   for(int j = 0; j < DSTTree->GetNtrig(); j++)
                   {
                        int tel_id = DSTTree->GetTrigI(j);
                        for(int p = 0; p < npix[tel_id]; p++)
                        {
                            pe_list[tel_id][p] = DSTTree->pe_list[tel_id][p];
                        }
                        
                   } 
                    
                }
            }
            else //if we have the fadc data, we can get pe from it
            {
                
            }
            if(do_image_clean)
            {
                image_clean(pe_list, pixel_neighbors, tail_cuts, pixel_in_image, DSTTree->GetNtrig(), DSTTree->GetTrig_List(), npix);
            }
            compute_moments(image, pe_list, pixel_in_image, ntel, DSTTree->GetTrig_List(), x_pix, y_pix);
            image->ConvertToRad(focal_length);
            rec->RecShower(image, cuts_option);

            
        }

    }





}