#ifndef _Detect_
#define _Detect_
#include "TTree.h"
#include "Limits_defined.h"
#include "io_hess.h"

class Detect
{
    public:
        TTree* detect_tree;
        int  Ntel;
        int Tel_id;
        float Pos[3];
        float Focal_length;
        float Effective_focal_length;
        int Npix;
        int Npix_disabled;
        int Ngains;
        float X_PixMM[LACT_MAXPIXELS];
        float Y_PixMM[LACT_MAXPIXELS];
        float R_PixMM[LACT_MAXPIXELS];
        int Pixel_Shape[LACT_MAXPIXELS];
        float Pixel_Size[LACT_MAXPIXELS];
        int NMirrors;
        float Mirror_Area;

        Detect();
        ~Detect();
        void InitWrite();
        void InitRead(TTree* t);
        void Reset();
        TTree* Get_DetectTree()
        {
            return detect_tree;
        }
        int GetNtel()
        {
            return Ntel;
        }
        int GetTel_id()
        {
            return Tel_id;
        }
        int GetNpix()
        {
            return Npix;
        }
        float GetFocal_Length()
        {
            return Focal_length;
        }
        float GetEffective_Focal_length()
        {
            return Effective_focal_length;
        }
        bool Fill_Data(AllHessData *hsdata);


};

















#endif