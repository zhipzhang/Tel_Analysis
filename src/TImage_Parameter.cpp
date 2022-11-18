#include "TImage_Parameter.h"
#include "rec_tools.h"
#include "TMath.h"
ClassImp(TImage_Parameter)

TImage_Parameter::TImage_Parameter(/* args */)
{
    for(int i = 0; i < LACT_MAXTEL; i++)
    {
        rp[i] = width[i] = length[i] = alpha[i] = image_x[i] = image_y[i] = 0.;
    }
    MRSW = MRSL = 0.;
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
        focal_length[itel] = f_length[itel];
    }
}

void TImage_Parameter::ComputeTelRp(double tel_pos[][3], TMcData* mc)
{
    for(int i = 0; i < image_tel.size(); i++)
    {
        int itel = image_tel[i];
        rp[itel] = line_point_distance(mc->core_pos[0], mc->core_pos[1], 0, cos(mc->true_direction[0] * TMath::DegToRad()) * cos(mc->true_direction[1] * TMath::DegToRad()),
                                        -sin(mc->true_direction[0] * TMath::DegToRad()) * cos(mc->true_direction[1] * TMath::DegToRad()), sin(mc->true_direction[1] * TMath::DegToRad() )
                                        ,tel_pos[itel][0], tel_pos[itel][1], tel_pos[itel][2]);
        
    }

}

void TImage_Parameter::clear()
{
    for(int i = 0; i < LACT_MAXTEL; i++)
    {
        rp[i] = width[i] = length[i] = alpha[i] = image_x[i] = image_y[i] = 0.;
    }
    MRSW = MRSL = 0.;
    image_tel.clear();

}