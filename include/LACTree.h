#ifndef LACTree_H
#define LACTree_H

#include "TH1F.h"
#include "TMath.h"
#include "TTree.h"

#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdint.h>
#include <string>
#include <vector>

#include "Limits_defined.h"

using namespace std;

class LACTree
{
    public:

        TTree* Mctree;
        TTree* EventTree;

    //  Mc_tree prameters
        unsigned short int primary;
        float energy;
        float xcore;
        float ycore;
        float az;
        float ze;
        float altitude;
        float Tel_az;
        float Tel_ze;   
        float Tel_al;
        unsigned int runNumber;
        unsigned int eventNumber; 
        float Point_Az[LACT_MAXTEL];
        float Point_Al[LACT_MAXTEL];
        
        //image parameter not used by now
        
        float width[LACT_MAXTEL];
        float length[LACT_MAXTEL];
        float x_img[LACT_MAXTEL];
        float y_img[LACT_MAXTEL];
        
        bool fMC;

        float weight; //weight for the different power-index  
        bool flag;
        unsigned int ntel;
        unsigned int LTrig;
        unsigned int Ntrig;
        unsigned int LTrig_list[LACT_MAXTEL];
        float LTime[LACT_MAXTEL];
        unsigned int ntel_data; //number of tel expected to have data
        unsigned int tel_data[LACT_MAXTEL];
        float size[LACT_MAXTEL]; //total pe
        float Rp[LACT_MAXTEL];
        float dsp_error[LACT_MAXTEL];
        float phi[LACT_MAXTEL]; // the angle of major axis

        //ADC parameters
        float fadc_pedestal[LACT_MAXTEL][LACT_MAXPIXELS];
        float fadc_sums[LACT_MAXTEL][LACT_MAXPIXELS];
        unsigned short int fadc_sum_window[LACT_MAXTEL][LACT_MAXPIXELS];
        unsigned short int fadc_sum_first[LACT_MAXTEL][LACT_MAXPIXELS];

        float pulse_timinglevels[LACT_MAXTEL][LACT_MAXPIXELS];
        float pulse_timing[LACT_MAXTEL][LACT_MAX_TIMELEVELS][LACT_MAXPIXELS];
        short int raw_max[LACT_MAXTEL][LACT_MAXPIXELS];

        bool fadc_read_write;   //whether write the samples of the fadc
        unsigned short int fadc_num_samples[LACT_MAXTEL];
        unsigned short int fadc_trace[LACT_MAXTEL][LACT_SUMWINDOW][LACT_MAXPIXELS];
        unsigned short int fadc_HG[LACT_MAXTEL][LACT_MAXPIXELS];
        unsigned short int Telescope_ZeroSuppression[LACT_MAXTEL];

        bool fillPeLeaf;
        float pe_list[LACT_MAXTEL][LACT_MAXPIXELS];




    // 
        LACTree();
        ~LACTree() {}

    // init Mctree (get Branch)  
        bool initMctree();
        TTree* getMCTree()
        {
            return Mctree;
        }
        TTree* getEventTree()
        {
            return EventTree;
        }
        bool isMC()
        {
            return fMC;
        }


        bool initEventTree(bool Full_tree, bool Calibration_tree);
        bool initEventTree(TTree* t);
        void resetDataVectors(unsigned int NMAX_LACT_MAXTEL=LACT_MAXTEL,unsigned int NMAX_LACT_MAXPIXELS=LACT_MAXPIXELS);

        void setFADC(bool iFADC = true)
        {
            fadc_read_write = iFADC;
        }
        void setFillPEleaf( bool iPE = false)
        {
            fillPeLeaf = iPE;
        }
        void setMC(bool iMC = true)
        {
            fMC = iMC;
        }
        unsigned int GetNtrig()
        {
            return Ntrig;
        }
        unsigned int GetTrigI(int i)
        {
            if( i < Ntrig)
            {
                return LTrig_list[i];
            }
            else
            {
                return -1;
            }

        }
        unsigned int* GetTrig_List()
        {
            return LTrig_list;
        }
};











#endif
