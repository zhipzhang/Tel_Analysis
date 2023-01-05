
#ifndef _THISTS_
#define _THISTS_

//#include "TObject.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TProfile.h"
#include <string>

class THists 
{
    public:
        TH2D  *h100[11]; // core position distribution in different energy range
        TH1D  *h1, *h2, *h3, *h5, *h6; //Different Number of Events : All, Trigger , Image > 4
        TH2D  *h4; // Array core distribution versus E image > 4
        TProfile  *h200[11]; // dist versus Rp in Different energy range
        TProfile  *h250[11];
        TProfile  *h2250[11];
        TProfile  *h2000[11]; // MISS Versus Rp in Different energy range
        TH2D  *h301, *h302; // MRSW:MRSL all / angular < 1 deg
        TH1D  *h303, *h304; // Total distribution of MRSW, MRSL
        TH1D  *h400[11]; // MRSW distribution in energy bins;
        TProfile  *h450[11]; // Width Distribution in energy bins;
        TH1D  *h500[11]; //MRSL distribution in energy bins;
        TProfile  *h550[11]; // Length Distribution in Energy bins;
      //  TH1D  *h600[10]; //MRSW:MRSW for different ntel
        TH1D  *h60; // Direction Error distribution
        TH1D  *h61; // Direction Error distribution (Pass Shape_cut)
        TProfile  *h62;
        TProfile  *h63;
        TH1D  *h64;
        TProfile  *h600[11]; //Direction Error versus Core distance in energy bins
        TH3D  *hist_max;
        TFile * histfile;
        THists(std::string& );

        ~THists();
        void Write();

};
















#endif