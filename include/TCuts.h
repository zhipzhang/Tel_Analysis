//Choose the cuts option when do the Shower Reconstrution
#ifndef _T_CUTS_
#define _T_CUTS_
class TCuts
{
    public:
        float min_size;
        float max_rp;
        float max_dist; // the core in nomial plane
        float energy_range[2];

        TCuts();
        ~TCuts();
        void SetMinSize(float min_amp)
        {
            min_size = min_amp;
        }
};

#endif