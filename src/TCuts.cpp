#include "TUserCuts.h"
#include "TMath.h"


TUserCuts::TUserCuts()
{
    min_size = 200.;
    max_rp = 9999.;
    max_dist = 4.;
    energy_range[0] = 0.;
    energy_range[1] = 1e9;
}

bool TUserCuts::Check(TImage_Parameter* image)
{
    int n = image->image_tel.size();
    if( n < min_tel)
    {
        return 0;
    }
    int pass = 0;
    for ( int i = 0; i < n ; i++)
    {
        int tel_id = image->image_tel[i];
        if(image->GetTelSize(tel_id) > min_size)
        {
            pass++;
            if(image->GetDist(tel_id) * TMath::RadToDeg() < max_dist )
            {
                image->SetSigTel(tel_id);
            }
        }
        
    }
    if(pass < min_tel || image->nsig < 2)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}