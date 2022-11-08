// MC parameters

#ifndef MonteCarloRunHeader_H
#define MonteCarloRunHeader_H

#include <TNamed.h>

#include <cmath>
#include <bitset>
#include <iostream>

using namespace std;

class MonteCarloRunHeader : public TNamed
{
    public:
    
        long int runnumber;
        int shower_date;         ///< date of shower simulations
        int detector_date;       ///< date of detector simulations
        /* ... shower MC   */
        unsigned int primary_id; ///< corsika ID
        double obsheight;        ///< Height of simulated observation level.
        int num_showers;         ///< Number of showers simulated.
        int num_use;             ///< Number of uses of each shower.
        int core_pos_mode;       ///< Core position fixed/circular/rectangular/...
        double core_range[2];    ///< rmin+rmax or dx+dy.
        double az_range[2];      ///< Range of shower azimuth [rad, N->E].
        double alt_range[2];     ///< Range of shower altitude [rad].
        double E_range[2];       ///< Energy range [TeV] of simulated showers.
        double spectral_index;   ///< Power-law spectral index of spectrum (<0).
        double B_total;          ///< Total geomagnetic field assumed [microT].
        double B_inclination;    ///< Inclination of geomagnetic field [rad].
        double B_declination;    ///< Declination of geomagnetic field [rad].
        double injection_height; ///< Height of particle injection [m].
        double fixed_int_depth;  ///< Fixed depth of first interaction or 0 [g/cm^2].
        int atmosphere;          ///< Atmospheric model number.
        /* ... + shower MC specific ... */
        double corsika_bunchsize;
        double corsika_wlen_min;
        double corsika_wlen_max;
        int corsika_low_E_detail;
        int corsika_high_E_detail;
        double corsika_low_high_E;    // transition energy
        /* ... + detector MC specific ... */
        string detector_Simulator;
        bool combined_runHeader;  // incomplete run header from several simtel input files
        
        float fFADC_hilo_multipler;
        
        MonteCarloRunHeader();
        ~MonteCarloRunHeader() {}
        double getMeanZenithAngle_Deg();
        void   print();
        void   printMCAz( bool iLowerLimit = true );
        void   printRunNumber();
        void   reset();
        bool   VOLUMEDET_set();
        
        ClassDef( MonteCarloRunHeader, 9 );
};
#endif
