/** CTA.convert_hessio_to_VDST
 *  short a program to convert sim_telarray files (hessio) to EVNDISP DST format
 *
 *
 *  Author of skeleton (as part of the hessio distribution):  Konrad Bernloehr
 *
 *  Author of modifications for eventdisplay: Gernot Maier (DESY)
 */


#include "initial.h"
#include "io_basic.h"
#include "history.h"
#include "io_hess.h"
#include "fileopen.h"
#include "rec_tools.h"

#include <bitset>
#include <iostream>
#include <map>
#include <set>
#include <stdlib.h>
#include <string>
#include <vector>
#include <signal.h>

#include "TFile.h"
#include "TList.h"
#include "TMath.h"
#include "TTree.h"
#include <TStopwatch.h>
#include "MonteCarloRunHeader.h"

#include "LACTree.h"
#include "Detect_config.h"
#include "TCalibData.h"
#include "LACTTeldata.h"
#include "LACTEvent.h"
#include "HitPix.h"

bool fFillPELeaf = false;
bool Only_High_Gain = true;
bool NO_fadc = false;

static int interrupted;
void stop_signal_function(int isig);
void stop_signal_function( int isig )
{
    if( isig >= 0 )
    {
        fprintf( stderr, "Received signal %d\n", isig );
    }
    if( !interrupted )
        fprintf( stderr,
                 "Program stop signal. Stand by until current data block is finished.\n" );
                 
    interrupted = 1;
    
    signal( SIGINT, SIG_DFL );
    signal( SIGTERM, SIG_DFL );
}
void syntax( char* program );

void syntax( char* program)
{
    printf("Program %s Usage : \n", program);
    printf("--root_file choose the ROOT output file (default: dst.root) \n");
    printf("--out_file  choose the Eventio output file (not used right now) \n");
    printf("--power-low power law index (will used when weighted) \n");
    printf("--no_fadc  if it's true , no fadc data will be stored! \n");
    printf("--write_pe stored the photoelectrons hit at each pixel \n");
}



/* fill MC run header*/
bool DST_fillMCRunheader(MonteCarloRunHeader*f, AllHessData* hsdata,bool multifiles);

bool DST_fillMCRunheader(MonteCarloRunHeader*f, AllHessData* hsdata,bool multifiles)
{   if(!multifiles)
    {
        f->runnumber = hsdata->run_header.run;
        f->obsheight = hsdata->mc_run_header.obsheight;
    }
    if(multifiles)
    {
        f->num_showers += hsdata->mc_run_header.num_showers;
    }
    else
    {
        f->num_showers = hsdata->mc_run_header.num_showers;
    }
    f->num_use = hsdata->mc_run_header.num_use;
    f->core_pos_mode = hsdata->mc_run_header.core_pos_mode;
    f->core_range[0] = hsdata->mc_run_header.core_range[0];
    f->core_range[1] = hsdata->mc_run_header.core_range[1];
    f->az_range[0] = hsdata->mc_run_header.az_range[0];
    f->az_range[1] = hsdata->mc_run_header.az_range[1];
    f->E_range[0] = hsdata->mc_run_header.E_range[0];
    f->E_range[1] = hsdata->mc_run_header.E_range[1];
    f->spectral_index = hsdata->mc_run_header.spectral_index;
    f->corsika_bunchsize = hsdata->mc_run_header.corsika_bunchsize;
    f->corsika_wlen_max = hsdata->mc_run_header.corsika_wlen_max;
    f->corsika_wlen_min = hsdata->mc_run_header.corsika_wlen_min;
    return true;
}
/* fill MC event Data (primary type, MC energy)*/
bool DST_fillMCEvent(LACTree* fData, AllHessData* hsdata);

bool DST_fillMCEvent(LACTree* fData, AllHessData* hsdata)
{
    fData->eventNumber = hsdata->mc_event.event;
    fData->runNumber = hsdata->run_header.run;
    fData->primary = hsdata->mc_shower.primary_id;
    fData->energy = hsdata->mc_shower.energy;
    fData->az = hsdata->mc_shower.azimuth * TMath::RadToDeg();
    fData->ze = 90. - hsdata->mc_shower.altitude * TMath::RadToDeg(); 
    fData->xcore = hsdata->mc_event.xcore;
    fData->ycore = hsdata->mc_event.ycore;
    if(fData->Mctree)
    {
        fData->Mctree->Fill();
    }
    return true;
    
}


bool DST_fillEvent(LACTree* fData, AllHessData* hsdata, double plidx);

bool DST_fillEvent(LACTree* fData, AllHessData* hsdata, double plidx)
{
    if(!fData || !hsdata)
    {
        return false;
    }

    fData->ntel = (unsigned short int)hsdata->run_header.ntel;
    fData->resetDataVectors(LACT_MAXTEL,LACT_MAXPIXELS);
    fData->eventNumber = hsdata->mc_event.event;
    fData->runNumber = hsdata->run_header.run;
 
    //only write MC data frist 
    fData->primary = hsdata->mc_shower.primary_id;
    fData->energy = hsdata->mc_shower.energy;
    fData->az = hsdata->mc_shower.azimuth * TMath::RadToDeg();
    fData->ze = 90. - hsdata->mc_shower.altitude * TMath::RadToDeg(); 
    fData->altitude = hsdata->mc_shower.altitude * TMath::RadToDeg();
    fData->xcore = hsdata->mc_event.xcore;
    fData->ycore = hsdata->mc_event.ycore;  
    double lg_E = log10(fData->energy);
    double index = hsdata->mc_run_header.spectral_index;
    fData->weight = pow(fData->energy, (plidx - index));

    //check flag whether it is trigger;
    //trigger data

    fData->Ntrig = hsdata->event.central.num_teltrg;
    unsigned int i_ntel_trig = 0;
    bitset<8 * sizeof(unsigned long)> i_localTrigger;//use bit to represent the trigger_list?

    //loop over all trigger event(not write now)
    for(unsigned int t=0; t<(unsigned int )hsdata->event.central.num_teltrg; t++)
    {
        if(hsdata->event.central.teltrg_list[t] < (int)i_localTrigger.size())   
        {
            i_localTrigger.set(hsdata->event.central.teltrg_list[t] - 1, true);
        }
        if(t < (unsigned int)hsdata->run_header.ntel)
        {
            fData->LTrig_list[i_ntel_trig] = hsdata->event.central.teltrg_list[t];
            fData->LTime[i_ntel_trig] = hsdata->event.central.teltrg_time[t];
            fData->length[i_ntel_trig] = hsdata->event.teldata[t].img->l;
            fData->width[i_ntel_trig] = hsdata->event.teldata[t].img->w;
            fData->x_img[i_ntel_trig] = hsdata->event.teldata[t].img->x;
            fData->y_img[i_ntel_trig] = hsdata->event.teldata[t].img->y;
            fData->size[i_ntel_trig] = hsdata->event.teldata[t].img->amplitude;
        }
        i_ntel_trig++;
    }

    // tracking data
    unsigned int z =0;
    for( unsigned int i=0; i<(unsigned short int)hsdata->run_header.ntel;i++)
    {
        fData->Point_Az[z] = hsdata->event.trackdata[i].azimuth_raw * TMath::RadToDeg();
        fData->Point_Al[z] = hsdata->event.trackdata[i].altitude_raw * TMath::RadToDeg();
        z++;
    }

    //event data 
    fData->ntel_data = (unsigned short int)hsdata->event.central.num_teldata;
    unsigned short int i_ntel_data = 0;

    assert(fData->ntel_data > 0 );
    //loop over all telescopes
    for( unsigned short int i = 0; i<fData->ntel_data; i++)
    {
        fData->tel_data[i_ntel_data] = (unsigned short int)hsdata->event.central.teldata_list[i];
        unsigned int telID = fData->tel_data[i_ntel_data]-1;
        
        fData->fadc_num_samples[i_ntel_data]  = (unsigned short int)hsdata->event.teldata[telID].raw->num_samples;
        fData->Telescope_ZeroSuppression[i_ntel_data] = (unsigned short int)hsdata->event.teldata[telID].raw->zero_sup_mode;

        if( hsdata->camera_set[telID].num_pixels >= LACT_MAXPIXELS)
        {
            cout << "Error: There are too much pixels in Tel" << telID << endl;
            cout << "The limit is "<< LACT_MAXTEL << " but there are " << hsdata->camera_set[telID].num_pixels << endl;  
        }
        // loop all pixels
        for(int p=0; p<hsdata->camera_set[telID].num_pixels; p++)
        {
            if(fFillPELeaf)
            {
                fData->pe_list[i_ntel_data][p] = hsdata->mc_event.mc_pe_list[telID].pe_count[p];
            }
            if(!NO_fadc)
            {
                fData->fadc_HG[i_ntel_data][p] = (unsigned short int) (hsdata->event.teldata[telID].raw->adc_sum[HI_GAIN][p]);

                if( hsdata->event.teldata[telID].raw->adc_sample && hsdata->event.teldata[telID].raw->adc_sample[HI_GAIN] 
                        && hsdata->event.teldata[telID].raw->num_samples)
                {
                    if( hsdata->event.teldata[telID].raw->num_gains == 2)
                    {
                        Only_High_Gain = false;
                    }

                    // Now only treat the High Gain Situation , May Add Low_Gain next time
                        for( int  t = 0; t < hsdata->event.teldata[telID].raw->num_samples; t++)
                        {
                            fData->fadc_trace[i_ntel_data][t][p] = hsdata->event.teldata[telID].raw->adc_sample[HI_GAIN][p][t];
                        }

                }
                else // no trace data
                {
                    fData->fadc_sums[i_ntel_data][p] = hsdata->event.teldata[telID].raw->adc_sum[HI_GAIN][p] - hsdata->tel_moni[telID].pedestal[HI_GAIN][p];
                }

                //record pedestal 
                if(hsdata->tel_moni[telID].num_ped_slices > 0)
                {
                    fData->fadc_pedestal[i_ntel_data][p] = hsdata->tel_moni->pedestal[HI_GAIN][p]/(double) (hsdata->tel_moni[telID].num_ped_slices);
                }
            }
        }
        
        /*
        fData->Rp[i_ntel_data] = line_point_distance(fData->xcore, fData->ycore, 0.,cos(hsdata->mc_shower.altitude)*
        cos(hsdata->mc_shower.azimuth), cos(hsdata->mc_shower.altitude) * sin(-1 * hsdata->mc_shower.azimuth), sin(hsdata->mc_shower.altitude),
        hsdata->run_header.tel_pos[telID][0], hsdata->run_header.tel_pos[telID][1], hsdata->run_header.tel_pos[telID][2]);
        double theta = compute_true_dsp(hsdata->mc_shower.azimuth, hsdata->mc_shower.altitude, 
        fData->xcore, fData->ycore,hsdata->run_header.tel_pos[telID][0], hsdata->run_header.tel_pos[telID][1],
        hsdata->event.trackdata[telID].azimuth_raw, hsdata->event.trackdata[telID].altitude_raw);
        fData->dsp_error[i_ntel_data] = compute_axis_diff(theta, hsdata->event.teldata[telID].img[0].phi);
        */
        i_ntel_data++;
    }    


    if(fData->EventTree)
    {
        fData->EventTree->Fill();
    }
    fData->resetDataVectors();
    return true;

}
void FillEvent(LACTEvent*,  AllHessData*);
void FillEvent(LACTEvent* lactevent,  AllHessData* hsdata)
{
    lactevent->SetRunnumber(hsdata->run_header.run);
    lactevent->SetEventnumber(hsdata->mc_event.event);
    lactevent->SetPrimary(hsdata->mc_shower.primary_id);
    lactevent->SetEnergy(hsdata->mc_shower.energy);
    lactevent->SetAzimuth(hsdata->mc_shower.azimuth * TMath::RadToDeg());
    lactevent->SetAltitude(hsdata->mc_shower.altitude * TMath::RadToDeg());
    lactevent->SetCorex(hsdata->mc_event.xcore);
    lactevent->SetCorey(hsdata->mc_event.ycore);
    lactevent->SetRefaz(hsdata->run_header.direction[0] * TMath::RadToDeg());
    lactevent->SetRefal(hsdata->run_header.direction[1] * TMath::RadToDeg());
    for( unsigned int i = 0; i < (unsigned int)hsdata->event.central.num_teltrg; i++)
    {
        int tel_id = hsdata->event.central.teldata_list[i];
        lactevent->GetiTel(tel_id).SetTelid(tel_id);
        lactevent->GetiTel(tel_id).SetTelaz(hsdata->event.trackdata[tel_id].azimuth_raw * TMath::RadToDeg());
        lactevent->GetiTel(tel_id).SetTelze(hsdata->event.trackdata[tel_id].altitude_raw * TMath::RadToDeg());
        for(int p = 0; p < hsdata->camera_set[tel_id].num_pixels; p++)
        {
            double pe = hsdata->mc_event.mc_pe_list[tel_id].pe_count[p];
            if( pe > 0)
               lactevent->GetiTel(tel_id).AddHitpix(p, 0., pe);
            
        }

    }
}
using namespace std;
/*
    write the main program
*/

int main(int argc,char** argv)
{
    IO_BUFFER* iobuf =NULL;
    IO_ITEM_HEADER item_header;
    const char* input_fname=NULL;

    int itel, rc =0;
    int tel_id;
    int ntel_trg = 0, min_tel_trg = 0;
    int nev = 0, ntrg = 0;
    // power-law index for weight
    double plidx = -2.7;
    size_t events = 0, max_events = 0;
    int iarg;
    int dst_processing = 0;
    char* program = argv[0];
    string all_input_file;                //  write all input_name as a string
    int num_input;
    string dst_file = "dst.root";        // output dst file
    string iobuf_file=" ";                   // new eventio file 

    static AllHessData* hsdata;
     /* Catch INTerrupt and TERMinate signals to stop program */
    signal( SIGINT, stop_signal_function );
    signal( SIGTERM, stop_signal_function );

     
    if( ( iobuf = allocate_io_buffer( 1000000L ) ) == NULL )
    {
        cout << "Cannot allocate I/O buffer" << endl;
        exit( EXIT_FAILURE );
    }
    iobuf->max_length = 1000000000L;
    H_CHECK_MAX();


    while (argc >1)
    {
        if(strcmp(argv[1], "--root_file") == 0 && argc >2)
        {
            dst_file = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        else if(strcmp(argv[1], "-output_file") == 0 && argc > 2)
        {
            iobuf_file = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        else if(strcmp(argv[1], "-dst-process") == 0 )
        {
            dst_processing = 1;
            argc -= 1;
            argv += 1;
            continue;
        }
        else if(strcmp(argv[1], "-power-low") == 0 && argc >2)
        {
            plidx = atof(argv[2]);
            argc -= 2;
            argv += 2;
            continue;
        }
        else if(strcmp(argv[1], "--no_fadc") == 0 )
        {
            NO_fadc  = true;
            argc -= 1;
            argv += 1;
            continue;
        }
        else if(strcmp(argv[1], "--write-pe") == 0)
        {
            fFillPELeaf = true;
            argc -= 1;
            argv += 1;
            continue;
        }
        else
        {
            break;
        }

    }  

    if(NO_fadc)
    {
        std::cout << "Only handle the ideal PE data and ignore the electronics So you must can read the pe data" << std::endl;
        
    }
    if( NO_fadc && !fFillPELeaf )
    {
        std::cout << "You have to have one of the ADC data or Pure PE data" << std::endl;
        exit(EXIT_FAILURE);
    }
 
    // open a root file to store data
    TFile* root_file = new TFile( dst_file.c_str(), "RECREATE" );
    if( root_file->IsZombie())
    {
        std::cout << "Error while opening root file " << dst_file << std::endl;
        exit(1);
    }
    LACTEvent* lactevent = new LACTEvent(); 
    TTree* event_tree = new TTree("LACTevent", "events data");
    event_tree->Branch("event", &lactevent);

    // open the new event_io file which name is iobuf_name
    // Maybe Can store the histogram using the eventio format
   /* if(iobuf_file != NULL && iobuf->output_file == NULL)
    {
        iobuf->output_file = fileopen(iobuf_file.c_str(), "w");
        if (iobuf->output_file == NULL)
        {
            perror(iobuf_file.c_str());
            exit(1);
        }
        else
        {
            std::cout << "iobuf_file " << iobuf_file << "has been opened!" << std::endl;
        }
    }*/
    // new VDSTTree class
    LACTree* DST = new LACTree();
    DST->setMC();
    DST->setFADC(!NO_fadc); //if no_fadc is true , set FADC to False
    DST->setFillPEleaf(fFillPELeaf);
    DST->initMctree();
    DST->initEventTree();


    Detect *detect = new Detect();
    detect->InitWrite();
    // MC run header
    MonteCarloRunHeader* MC_header = new MonteCarloRunHeader();
    MC_header->SetName( "MC_runheader" );


    TCalibData* calib_data = new TCalibData();
    cout << "start proceeding with sim_telarray file" << endl;
    // Now go over rest of the command line and read input file names
    while( argc > 1 || input_fname != NULL )
    {
        if( interrupted )
        {
            break;
        }
        if( argc > 1 )
        {
            if( argv[1][0] == '-' && argv[1][1] != '\0' )
            {
                syntax( program );
                exit(EXIT_FAILURE);
                
            }
            else
            {
                input_fname = argv[1];
                argc--;
                argv++;
            }
        }
        
        if( ( iobuf->input_file = fileopen( input_fname, READ_BINARY ) ) == NULL )
        {
            perror( input_fname );
            cout << "Cannot open input file." << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        if( input_fname )
        {
            cout << "opening simtel file " << input_fname << endl;
            fflush( stdout );
            all_input_file += input_fname;
            num_input ++;
            if(num_input > 0)
            {
                all_input_file += ",  ";
            }
        }
        input_fname = NULL;
        
        for( ;; )
        {
            if( interrupted )
            {
                break;
            }
            
            /* Find and read the next block of data. */
            /* In case of problems with the data, just give up. */
            if( find_io_block( iobuf, &item_header ) != 0 )
            {
                break;
            }
            if( read_io_block( iobuf, &item_header ) != 0 )
            {
                break;
            }
            
            if( hsdata == NULL &&
                    item_header.type > IO_TYPE_HESS_RUNHEADER &&
                    item_header.type < IO_TYPE_HESS_RUNHEADER + 200 )
            {
                fprintf( stderr, "Trying to read event data before run header.\n" );
                fprintf( stderr, "Skipping this data block.\n" );
                continue;
            }
            
            //////////////////////////////////////////
            // check header types
            switch( ( int ) item_header.type )
            { /* =================================================== */
                case IO_TYPE_HESS_RUNHEADER:
                    /* Summary of a preceding run in the same file ? */
                    if( nev > 0 )
                    {
                        printf( "%d of %d events triggered.\n", ntrg, nev );
                    }
                    /* Structures might be allocated from previous run */
                    if( hsdata != NULL )
                    {
                        /* Free memory allocated inside ... */
                        for( itel = 0; itel < hsdata->run_header.ntel; itel++ )
                        {
                            if( hsdata->event.teldata[itel].raw != NULL )
                            {
                                free( hsdata->event.teldata[itel].raw );
                                hsdata->event.teldata[itel].raw = NULL;
                            }
                            if( hsdata->event.teldata[itel].pixtm != NULL )
                            {
                                free( hsdata->event.teldata[itel].pixtm );
                                hsdata->event.teldata[itel].pixtm = NULL;
                            }
                            if( hsdata->event.teldata[itel].img != NULL )
                            {
                                free( hsdata->event.teldata[itel].img );
                                hsdata->event.teldata[itel].img = NULL;
                            }
                        }
                        /* Free main structure */
                        free( hsdata );
                        hsdata = NULL;
                    }
                    hsdata = ( AllHessData* ) calloc( 1, sizeof( AllHessData ) );
                    if( ( rc = read_hess_runheader( iobuf, &hsdata->run_header ) ) < 0 )
                    {
                        cout << "Reading run header failed." << endl;
                        exit( EXIT_FAILURE );
                    }
                    cout << "Telescope configuration: " << endl;
                    for( itel = 0; itel < hsdata->run_header.ntel; itel++ )
                    {
                        tel_id = hsdata->run_header.tel_id[itel];
                        hsdata->camera_set[itel].tel_id = tel_id;
                        cout << "\t initialize Telescope ID " << tel_id << " (is telescope # " << itel << ")" << endl;
                        hsdata->camera_org[itel].tel_id = tel_id;
                        hsdata->pixel_set[itel].tel_id = tel_id;
                        hsdata->pixel_disabled[itel].tel_id = tel_id;
                        hsdata->cam_soft_set[itel].tel_id = tel_id;
                        hsdata->tracking_set[itel].tel_id = tel_id;
                        hsdata->point_cor[itel].tel_id = tel_id;
                        hsdata->event.num_tel = hsdata->run_header.ntel;
                        hsdata->event.teldata[itel].tel_id = tel_id;
                        hsdata->event.trackdata[itel].tel_id = tel_id;
                        if( ( hsdata->event.teldata[itel].raw =
                                    ( AdcData* ) calloc( 1, sizeof( AdcData ) ) ) == NULL )
                        {
                            cout << "Not enough memory" << endl;
                            exit( EXIT_FAILURE );
                        }
                        hsdata->event.teldata[itel].raw->tel_id = tel_id;
                        if( ( hsdata->event.teldata[itel].pixtm =
                                    ( PixelTiming* ) calloc( 1, sizeof( PixelTiming ) ) ) == NULL )
                        {
                            cout << "Not enough memory" << endl;
                            exit( EXIT_FAILURE );
                        }
                        hsdata->event.teldata[itel].pixtm->tel_id = tel_id;
                        if( ( hsdata->event.teldata[itel].img =
                                    ( ImgData* ) calloc( 2, sizeof( ImgData ) ) ) == NULL )
                        {
                            cout << "Not enough memory" << endl;
                            exit( EXIT_FAILURE );
                        }
                        hsdata->event.teldata[itel].max_image_sets = 2;
                        hsdata->event.teldata[itel].img[0].tel_id = tel_id;
                        hsdata->event.teldata[itel].img[1].tel_id = tel_id;
                        hsdata->tel_moni[itel].tel_id = tel_id;
                        hsdata->tel_lascal[itel].tel_id = tel_id;
                    }
                    break;
                    // all things do in runheader
                case IO_TYPE_HESS_MCRUNHEADER:
                    rc = read_hess_mcrunheader( iobuf, &hsdata->mc_run_header );
                    // fill EVNDISP DST run header
                    DST_fillMCRunheader( MC_header, hsdata ,(num_input > 1));
                    break;
                    
                /* =================================================== */
                case IO_TYPE_MC_INPUTCFG:
                {
                    struct linked_string corsika_inputs;
                    corsika_inputs.text = NULL;
                    corsika_inputs.next = NULL;
                    read_input_lines( iobuf, &corsika_inputs );
                    if( corsika_inputs.text != NULL )
                    {
                        struct linked_string* xl = NULL, *xln = NULL;
                        for( xl = &corsika_inputs; xl != NULL; xl = xln )
                        {
                            free( xl->text );
                            xl->text = NULL;
                            xln = xl->next;
                            xl->next = NULL;
                            if( xl != &corsika_inputs )
                            {
                                free( xl );
                            }
                        }
                    }
                }
                break;
                
                /* =================================================== */
                case 70: /* How sim_hessarray was run and how it was configured. */
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_CAMSETTINGS:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Camera settings for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_camsettings( iobuf, &hsdata->camera_set[itel] );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_CAMORGAN:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Camera organisation for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_camorgan( iobuf, &hsdata->camera_org[itel] );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_PIXELSET:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Pixel settings for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_pixelset( iobuf, &hsdata->pixel_set[itel] );
                    break; 
                /* =================================================== */
                case IO_TYPE_HESS_PIXELDISABLE:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Pixel disable block for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_pixeldis( iobuf, &hsdata->pixel_disabled[itel] );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_CAMSOFTSET:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Camera software settings for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_camsoftset( iobuf, &hsdata->cam_soft_set[itel] );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_POINTINGCOR:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Pointing correction for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_pointingcor( iobuf, &hsdata->point_cor[itel] );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_TRACKSET:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Tracking settings for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_trackset( iobuf, &hsdata->tracking_set[itel] );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_EVENT:
                    rc = read_hess_event( iobuf, &hsdata->event, -1 );
                    events++;
                    /* Count number of telescopes (still) present in data and triggered */
                    ntel_trg = 0;

                    for( itel = 0; itel < hsdata->run_header.ntel; itel++ )
                    {
                        if( hsdata->event.teldata[itel].known )
                        {
                            /* If non-triggered telescopes record data (like HEGRA),
                               we may have to check the central trigger bit as well,
                               but ignore this for now. */
                            ntel_trg++;
                        }
                    }
                    if( hsdata->event.shower.known )
                    {
                        hsdata->event.shower.num_trg = ntel_trg;
                    }
                    if( ntel_trg < min_tel_trg )
                    {
                        continue;
                    }
                    ntrg++;
                    
                    DST_fillEvent( DST, hsdata ,plidx);
                    FillEvent(lactevent,  hsdata);
                    event_tree->Fill();
                    lactevent->clear();
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_CALIBEVENT:
                {
                    int type = -1;
                    rc = read_hess_calib_event( iobuf, &hsdata->event, -1, &type );
                    /*if( !fWriteTelConfigTreeOnly )
                    {
                        DST_fillEvent( fDST, hsdata, fTelescope_list, fWriteFADC );
                    }*/
                }
                break;
                
                /* =================================================== */
                case IO_TYPE_HESS_MC_SHOWER:
                    rc = read_hess_mc_shower( iobuf, &hsdata->mc_shower );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_MC_EVENT:
                    rc = read_hess_mc_event( iobuf, &hsdata->mc_event );
                    DST_fillMCEvent( DST, hsdata);
                    
                    break;
                    
                /* =================================================== */
                case IO_TYPE_MC_TELARRAY:
                    if( hsdata && hsdata->run_header.ntel > 0 )
                    {
                        rc = read_hess_mc_phot( iobuf, &hsdata->mc_event );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_MC_PE_SUM:
                    rc = read_hess_mc_pe_sum( iobuf, &hsdata->mc_event.mc_pesum );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_TEL_MONI:
                    // Telescope ID among others in the header
                    tel_id = ( item_header.ident & 0xff ) |
                             ( ( item_header.ident & 0x3f000000 ) >> 16 );
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Telescope monitor block for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_tel_monitor( iobuf, &hsdata->tel_moni[itel] );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_LASCAL:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Laser/LED calibration for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_laser_calib( iobuf, &hsdata->tel_lascal[itel] );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_RUNSTAT:
                    rc = read_hess_run_stat( iobuf, &hsdata->run_stat );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_MC_RUNSTAT:
                    rc = read_hess_mc_run_stat( iobuf, &hsdata->mc_run_stat );
                    break;
                    
                default:
                        fprintf( stderr, "Ignoring data block type %ld\n", item_header.type );
                        break; 
                              
            }
        }

        if(iobuf->input_file!=NULL)
        {
            fileclose(iobuf->input_file);
        }
        iobuf->input_file = NULL;
        reset_io_block(iobuf);

        if(hsdata != NULL)
        {
            hsdata->run_header.run = 0;
        }
    }
    calib_data->InitWrite();
    calib_data->FillData(hsdata);
    if(calib_data->GetCalibTree())
    {
        calib_data->GetCalibTree()->Write();
    }

    if(DST && DST->getMCTree())
    {
        DST->getMCTree()->Write();
    }
    detect->Fill_Data(hsdata);
    if(DST->getEventTree())
    {
        DST->getEventTree()->Write();
    }
    if(detect->Get_DetectTree())
    {
        detect->Get_DetectTree()->Write();
    }
    event_tree->Write();

    if(root_file)
    {
        root_file->Close();
    }
    cout <<"written to DST file "<<endl;
    return 0;
}
