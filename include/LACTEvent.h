#ifndef _LACTEVENT_H_
#define _LACTEVENT_H_

#include "TObject.h"
#include "LACTTeldata.h"
#include "TObjArray.h"
#include "Limits_defined.h"
#include <vector>

class LACTEvent : public TObject
{
    private:
        Int_t runnumber;
        Int_t eventnumber;
        Int_t NTrigTel;
        Int_t primary_id;
        Double_t energy;
        Double_t azimuth;
        Double_t altitude;
        Double_t ref_az; //reference tel point azimuth
        Double_t ref_al; // the same as above
        Double_t corex;
        Double_t corey;

        LACTTeldata TelData[LACT_MAXTEL];
    
    public:
        LACTEvent();
        virtual ~LACTEvent();
        void SetRunnumber(Int_t rnumber)
        {
            runnumber = rnumber;
        }
        void SetEventnumber(Int_t enumber)
        {
            eventnumber = enumber;
        }
        void SetNtrigTel(Int_t trig_num)
        {
            NTrigTel = trig_num;
        }
        void SetEnergy(Double_t E)
        {
            energy = E;
        }
        void SetPrimary(Int_t primary)
        {
            primary_id = primary;
        }
        void SetAzimuth(Double_t theta)
        {
            azimuth = theta;
        }
        void SetAltitude(Double_t phi)
        {
            altitude = phi;
        }
        void SetRefaz(Double_t ref_theta)
        {
            ref_az = ref_theta;
        }
        void SetRefal(Double_t ref_phi)
        {
            ref_al = ref_phi;
        }
        void SetCorex(Double_t x)
        {
            corex = x;
        }
        void SetCorey(Double_t y)
        {
            corey = y;
        }
        LACTTeldata GetiTel(Int_t id)
        {
            return TelData[id];
        }
        void clear()
        {
            for(int i = 0; i < LACT_MAXTEL; i++)
            {
                TelData[i].SetNhitpix(0);
                TelData[i].clear();
            }
        }
        ClassDef(LACTEvent, 1)
};
























#endif