#include "Detect_config.h"

Detect::Detect()
{
    Ntel = 0;
    Tel_id = 0;
    Effective_focal_length =  Focal_length = 0.;
    Npix = Npix_disabled = 0;
    Ngains = 0;
    NMirrors = 0;
    Mirror_Area = 0.;
    std::fill(Pos, Pos+3, 0);

}
Detect::~Detect()
{
    free(X_PixMM);
    free(Y_PixMM);
    free(R_PixMM);
    free(Pixel_Size);
}

void Detect::Reset()
{
    Ntel = 0;
    Tel_id = 0;
    Effective_focal_length =  Focal_length = 0.;
    Npix = Npix_disabled = 0;
    Ngains = 0;
    NMirrors = 0;
    Mirror_Area = 0.;
    std::fill(Pos, Pos+3, 0);
    std::fill(X_PixMM, X_PixMM + LACT_MAXPIXELS, 0.);
    std::fill(Y_PixMM, Y_PixMM + LACT_MAXPIXELS, 0.);
    std::fill(Pixel_Shape, Pixel_Shape + LACT_MAXPIXELS, -1);
    std::fill(Pixel_Size, Pixel_Size + LACT_MAXPIXELS, 0.);

}
void Detect::InitWrite()
{
    X_PixMM = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    Y_PixMM = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    R_PixMM = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    Pixel_Size = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    detect_tree = new TTree("config_tree", "config_data");
    detect_tree->Branch("Ntel", &Ntel, "Ntel/I");
    detect_tree->Branch("Tel_id", &Tel_id, "Tel_id/I");
    detect_tree->Branch("Tel_pos", Pos, "Tel_pos[3]/F");
    detect_tree->Branch("Focal_length", &Focal_length, "Focal_length/F");
    detect_tree->Branch("Effective_focal_length", &Effective_focal_length, "Effective_focal_length/F");
    detect_tree->Branch("Npix", &Npix, "Npix/I");
    detect_tree->Branch("Npix_disabled", &Npix_disabled, "Npix_active/I");
    detect_tree->Branch("Ngains", &Ngains, "Ngains/I");
    detect_tree->Branch("X_Pix", X_PixMM, "X_Pix[Npix]/F");
    detect_tree->Branch("Y_Pix", Y_PixMM, "Y_Pix[Npix]/F");
   // detect_tree->Branch("R_Pix", R_PixMM, "Y_Pix[Npix]/F");
    detect_tree->Branch("Pixel_Shape", Pixel_Shape, "Pixel_Shape[Npix]/I");
    detect_tree->Branch("Pixel_Size", Pixel_Size, "Pixel_Size[Npix]/F");
    detect_tree->Branch("NMirrors", &NMirrors, "NMirrors/I");
    detect_tree->Branch("Mirror_Area", &Mirror_Area, "Mirror_Area/F");

}
void Detect::InitRead(TTree *t)
{
    X_PixMM = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    Y_PixMM = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    R_PixMM = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    Pixel_Size = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    detect_tree = t;
    detect_tree->SetBranchAddress("Ntel", &Ntel);
    detect_tree->SetBranchAddress("Tel_id", &Tel_id);
    detect_tree->SetBranchAddress("Tel_pos", Pos);
    detect_tree->SetBranchAddress("Focal_length", &Focal_length);
    detect_tree->SetBranchAddress("Effective_focal_length", &Effective_focal_length);
    detect_tree->SetBranchAddress("Npix", &Npix);
    detect_tree->SetBranchAddress("Npix_disabled", &Npix_disabled);
    detect_tree->SetBranchAddress("Ngains", &Ngains);
    detect_tree->SetBranchAddress("X_Pix", X_PixMM);
    detect_tree->SetBranchAddress("Y_Pix", Y_PixMM);
    //detect_tree->SetBranchAddress("R_Pix", R_PixMM);
    detect_tree->SetBranchAddress("Pixel_Shape", Pixel_Shape);
    detect_tree->SetBranchAddress("Pixel_Size", Pixel_Size);
    detect_tree->SetBranchAddress("NMirrors", &NMirrors);
    detect_tree->SetBranchAddress("Mirror_Area", &Mirror_Area);

}


bool Detect::Fill_Data(AllHessData* hsdata)
{
    if( !hsdata )
    {
        printf("Error When Filling Config Tree \n" );
        return 0;
    }
    Ntel = hsdata->run_header.ntel;
    for( int itel = 0; itel < hsdata->run_header.ntel; itel++)
    {
        Tel_id = hsdata->run_header.tel_id[itel];
        Pos[0] = hsdata->run_header.tel_pos[itel][0];
        Pos[1] = hsdata->run_header.tel_pos[itel][1];
        Pos[2] = hsdata->run_header.tel_pos[itel][2];
        Focal_length = hsdata->camera_set[itel].flen;
        if(hsdata->camera_set[itel].eff_flen > 0)
        {
            Effective_focal_length = hsdata->camera_set[itel].eff_flen;
        }

        NMirrors = hsdata->camera_set[itel].num_mirrors;
        Mirror_Area = hsdata->camera_set[itel].mirror_area;
        Npix = hsdata->camera_set[itel].num_pixels;
        Npix_disabled = hsdata->pixel_disabled[itel].num_HV_disabled;
        Ngains = hsdata->event.teldata[itel].raw->num_gains;

        if(Npix > LACT_MAXPIXELS)
        {
            printf("There are too many pixels in Tel %d \n ", Tel_id);
            exit(EXIT_FAILURE);
        } 
        else
        {
            for(int p = 0; p < Npix; p++)
            {
                X_PixMM[p] = hsdata->camera_set[itel].xpix[p] * 1e3;
                Y_PixMM[p] = hsdata->camera_set[itel].ypix[p] * 1e3;
                Pixel_Shape[p] = hsdata->camera_set[itel].pixel_shape[p];
                Pixel_Size[p] = hsdata->camera_set[itel].size[p] * 1e3;

            }
        }
        detect_tree->Fill();


    }

}