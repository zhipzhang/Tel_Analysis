#include "TCuts.h"


TCuts::TCuts()
{
    min_size = 50.;
    max_rp = 9999.;
    max_dist = 9999.;
    energy_range[0] = 0.;
    energy_range[1] = 1e9;
}