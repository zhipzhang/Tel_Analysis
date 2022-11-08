#include "TCuts.h"


TCuts::TCuts()
{
    min_size = 0.;
    max_rp = 9999.;
    max_dist = 0.;
    energy_range[0] = 0.;
    energy_range[1] = 1e9;
}