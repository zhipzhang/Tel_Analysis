#include "LACTTeldata.h"
#include "HitPix.h"
#include "Limits_defined.h"
ClassImp(LACTTeldata)
ClassImp(HitPix)

LACTTeldata::LACTTeldata()
{
    Nhitpix = 0;
    HitsPix = new TClonesArray("HitPix", LACT_MAXPIXELS);
}

LACTTeldata::~LACTTeldata()
{

}

void LACTTeldata::AddHitpix(Int_t hitid, Double_t hit_time, Double_t hit_pe)
{
    new((*HitsPix)[Nhitpix ++]) HitPix(hitid, hit_time, hit_pe);
}



