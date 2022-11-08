#ifndef _CALIBDATA_
#define _CALIBDATA_

#include "TTree.h"
#include "Limits_defined.h"
#include "io_hess.h"

class TCalibData 
{
    public:
        TTree* calibration_tree;
        int Tel_id;
        int Npix;
        int Num_sumwindow;
        float High_pedstal[LACT_MAXPIXELS];
        float Low_pedstal[LACT_MAXPIXELS];
        float High_pedstal_var[LACT_MAXPIXELS];
        float Low_pedstal_var[LACT_MAXPIXELS];
        float High_Convert[LACT_MAXPIXELS];
        float Low_Convert[LACT_MAXPIXELS];
        bool Low_gain;

        TCalibData();
        ~TCalibData();
        void InitWrite();
        void InitRead();
        void Reset();
        void FillData(AllHessData *hsdata);
        TTree* GetCalibTree()
        {
            return calibration_tree;
        }

};

































#endif