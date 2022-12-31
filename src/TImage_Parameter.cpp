#include "TImage_Parameter.h"
#include "rec_tools.h"
#include "TMath.h"
ClassImp(TImage_Parameter)

TImage_Parameter::TImage_Parameter(/* args */)
{
    for(int i = 0; i < LACT_MAXTEL; i++)
    {
       rec_rp[i] = rp[i] = width[i] = length[i] = alpha[i] = image_x[i] = image_y[i] = 0.;
    }
}

TImage_Parameter::~TImage_Parameter()
{
}
void TImage_Parameter::ConvertToRad(double* f_length)
{
    for(int i = 0; i < image_tel.size(); i++)
    {
        int itel = image_tel[i];
        width[itel] = width[itel] /f_length[itel]   ;
        length[itel] = length[itel] / f_length[itel];
        image_x[itel] = image_x[itel]/f_length[itel] ;
        image_y[itel] = image_y[itel]/f_length[itel] ;
        source_x[itel] = source_x[itel]/f_length[itel];
        source_y[itel] = source_y[itel]/f_length[itel];
        focal_length[itel] = f_length[itel];
    }
}


void TImage_Parameter::clear()
{
    for(int i = 0; i < LACT_MAXTEL; i++)
    {
        rp[i] = width[i] = length[i] = alpha[i] = image_x[i] = image_y[i] = 0.;
    }
    nsig = 0;
    image_tel.clear();

}