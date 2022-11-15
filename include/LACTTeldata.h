#ifndef _LACT_TELDATA_H_
#define _LACT_TELDATA_H_

#include "TObject.h"
#include "TClonesArray.h"

class LACTTeldata : public TObject
{
    private:
        Int_t Tel_id;
        Int_t Nhitpix;
        Double_t Tel_az;
        Double_t Tel_al;
        Double_t Tel_pos[3];
        Double_t size;

        TClonesArray* HitsPix;
    public:
        LACTTeldata();
        virtual ~LACTTeldata();
        void SetTelid(Int_t id)
        {
            Tel_id = id;
        }
        void SetNhitpix(Int_t Nhit)
        {
            Nhitpix = Nhit;
        }
        void SetTelsize(Double_t pesize)
        {
            size = pesize;
        }
        void SetTelaz(Double_t azimuth)
        {
            Tel_az = azimuth;
        }
        void SetTelze(Double_t altitude)
        {
            Tel_al = altitude;
        }
        void SetTelpos(Double_t x, Double_t y, Double_t z)
        {
            Tel_pos[0] = x, Tel_pos[1] = y, Tel_pos[2] = z;
        }
        Int_t GetTelid()
        {
            return Tel_id;
        }
        Int_t GetNhitpix()
        {
            return Nhitpix;
        }
        Double_t GetTelSize()
        {
            return size;
        }
        Double_t GetTelaz()
        {
            return Tel_az;
        }
        Double_t GetTelal()
        {
            return Tel_al;
        }
        Double_t* GetTelpos()
        {
            return Tel_pos;
        }
        TClonesArray* GetHitspix()
        {
            return HitsPix;
        }
        void clear()
        {
            HitsPix->Clear();
        }


        void AddHitpix(Int_t id, Double_t time, Double_t pe);
        ClassDef(LACTTeldata,1)
        

};





#endif