#ifndef _LACT_PIX_H_
#define _LACT_PIX_H_

#include "TObject.h"

class HitPix :public TObject
{
    private:
        Int_t id;
        Double_t time;
        Double_t pe;
        Int_t status;
    public:
        HitPix();
        HitPix(Int_t hit_id, Double_t hit_time, Double_t hit_pe)
        {
            id = hit_id, time = hit_time, pe = hit_pe, status = 1;
        }
        virtual ~HitPix() {};
        Int_t GetId()
        {
            return id;
        }
        Double_t GetTime()
        {
            return time;
        }
        Double_t GetPe()
        {
            return pe;
        }
        Double_t GetStatus()
        {
            return status;
        }
        void SetId(Int_t hit_id)
        {
            id = hit_id;
        }
        void SetTime(Double_t hit_time)
        {
            time = hit_time;
        }
        void SetPe(Double_t hit_pe)
        {
            pe = hit_pe;
        }
        void SetStatus(Int_t Status)
        {
            status = Status;
        }
        ClassDef(HitPix, 1);

};



















#endif