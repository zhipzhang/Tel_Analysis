//Choose the cuts option when do the Shower Reconstrution
#ifndef _T_CUTS_
#define _T_CUTS_
#include "TImage_Parameter.h"
class TUserCuts
{
    public:
        float min_size;
        float max_rp;
        float max_dist; // the core in nomial plane
        float energy_range[2];
        int min_tel{4};
        float width_cut[2];
        float length_cut[2];

        TUserCuts();
        ~TUserCuts();
        void SetMinSize(float min_amp)
        {
            min_size = min_amp;
        }
        void SetWidthCut(float* w_cut)
        {
            width_cut[0] = w_cut[0];
            width_cut[1] = w_cut[1];
        }
        void SetLengthCut(float* l_cut)
        {
            length_cut[0] = l_cut[0];
            length_cut[1] = l_cut[1];
        }
        bool Check( TImage_Parameter* image);
};

#endif