// Main Class For Writing the analysis Tree
#ifndef _IMAGE_PARAMETER_
#define _IMAGE_PARAMETER_
#include "Limits_defined.h"
#include "TObject.h"
#include <vector>
#include "TMath.h"
class TImage_Parameter : public TObject
{
public: 
    /* data */
    std::vector<int> image_tel; // show which tel is in image very_important
    double rp[LACT_MAXTEL];
    double rec_rp[LACT_MAXTEL];
    double width[LACT_MAXTEL];
    double length[LACT_MAXTEL];
    double alpha[LACT_MAXTEL];
    double image_x[LACT_MAXTEL];  // All are stored in Rad;
    double image_y[LACT_MAXTEL];
    double source_x[LACT_MAXTEL];
    double source_y[LACT_MAXTEL];
    double size[LACT_MAXTEL];
    double focal_length[LACT_MAXTEL];
    bool have_lookup = false;
    int nsig = 0;                    // significant Tel : tel pass the cuts_option
    int sig_tel[LACT_MAXTEL];

    TImage_Parameter(/* args */);
    virtual ~TImage_Parameter();
    void Setlookup()
    {
        have_lookup = 1;
    }
    void SetTelWidth(int itel, double iwidth)
    {
        width[itel] = iwidth;
    }
    double GetTelwidth(int itel)
    {
        return width[itel];
    }
    void SetTelLength(int itel, double ilength)
    {
        length[itel] = ilength;
    }
    double GetTelLength(int itel)
    {
        return length[itel];
    }
    void SetTelAlpha(int itel, double ialpha)
    {
        alpha[itel] = ialpha;
    }
    double GetTelAlpha(int itel)
    {
        return alpha[itel];
    }
    void SetTelImageX(int itel, double imagex)
    {
        image_x[itel] = imagex;
    }
    double GetTelImageX(int itel)
    {
        return image_x[itel];
    }
    void SetTelImageY(int itel, double imagey)
    {
        image_y[itel] = imagey;
    }
    double GetTelImageY(int itel)
    {
        return image_y[itel];
    }
    double GetDist(int itel)
    {
        return sqrt(pow(image_x[itel], 2) + pow(image_y[itel] ,2));
    }
    void SetTelSize(int itel , double isize)
    {
        size[itel] = isize;
    }
    double GetTelSize(int itel)
    {
        return size[itel];
    }
    void SetTelRp(int itel, double irp)
    {
        rp[itel] = irp;
    }
    double GetTelRp(int itel)
    {
        return rp[itel];
    }
    double GetRecRp(int itel)
    {
        return rec_rp[itel];
    }
    double GetTelFocal(int itel)
    {
        return focal_length[itel];
    }
    void SetSigTel(int itel)
    {
        sig_tel[nsig++] = itel;
    }
    void SetTelSource(int itel, double sx, double sy)
    {
        source_x[itel] = sx;
        source_y[itel] = sy;
    }
    

    void ConvertToRad(double* focal_length);
    void clear();
    ClassDef(TImage_Parameter, 1)
};



#endif