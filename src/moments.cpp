#include "moments.h"
#include <cmath>
#include "TMath.h"

void image_clean(double pe_list[][LACT_MAXPIXELS], std::map<int, std::vector<int> > *pixel_neighbors, double* tail_cuts, 
                                                std::vector<int> *pixel_in_image, unsigned int ntel, unsigned int* Trig_List, int* npix)
{
    unsigned int tel_id;
    for(int i = 0; i < ntel; i++)
    {
        tel_id = Trig_List[i];
        for(int p = 0; p < npix[tel_id]; p++)
        {
            if( pe_list[tel_id][p] > tail_cuts[1])
            {
                for(auto k : pixel_neighbors[tel_id][p])
                {
                    if(pe_list[tel_id][k] > tail_cuts[0])
                    {
                        pixel_in_image[tel_id].push_back(p);
                    }
                    else
                    {
                        continue;
                    }
                }


            }
            else if(pe_list[tel_id][p] > tail_cuts[0])
            {
                for(auto k : pixel_neighbors[tel_id][p])
                {
                    if(pe_list[tel_id][k] > tail_cuts[1])
                    {
                        pixel_in_image[tel_id].push_back(p);
                    }
                }
            }
            else
            {
                continue;
            }
        }
    }
}


/*be careful alpha not consider the camera rotation !*/

void compute_moments(TImage_Parameter* image,double pe_list[][LACT_MAXPIXELS], std::vector<int>* pixel_in_image,unsigned int ntel, unsigned int* Trig_List,
                     double x_pix[][LACT_MAXPIXELS], double y_pix[][LACT_MAXPIXELS])
{
    int tel_id;
    for(int i = 0; i < ntel; i++)
    {
        double sx = 0, sxx = 0, sxy = 0, sy = 0, syy  = 0, sA = 0.;
        double a = 0, b = 0;
        tel_id = Trig_List[i];
        if(pixel_in_image[tel_id].size() < 3)
        {
            continue;
        }
        else
        {
            image->image_tel.push_back(tel_id);
        }
        for(int p = 0; p < pixel_in_image[tel_id].size(); p++)
        {
            int ipix = pixel_in_image[tel_id][p];
            double x = x_pix[tel_id][ipix];
            double y = y_pix[tel_id][ipix];
            double A = pe_list[tel_id][ipix];
            sA  += A;
            sx  += (A * x);
            sxx += (A * x) * x;
            sxy += (A * x) * y;
            sy  += (A * y);
            syy += (A * y) * y;

        }
        sx /= sA;
        sy /= sA;
        sxx = sxx/sA - sx * sx;
        sxy = sxy/sA - sx * sy;
        syy = syy/sA - sy * sy;
        if (fabs(sxy) > 1e-8 * fabs(sxx) && fabs(sxy) > 1e-8 * fabs(syy) )
        {
            double p1 = syy - sxx, p2 = sxy*sxy;
            double q, r1, r2;
            if ( p2 > 1e-8*(p1*p1) )
                q = p1 + sqrt(p1*p1+4.*p2);
            else
                q = 2.*p2;
            b = 0.5 * q/sxy; // solve  the equation
            a = sy - b*sx;
            if ( (r1 = syy + 2.*p2/q) > 0. )
                image->length[tel_id] = sqrt(r1);
            else
                image->length[tel_id] = 0.;
            if ( (r2 = sxx - 2.*p2/q) > 0. )
                image->width[tel_id] =  sqrt(r2);
            else
                image->width[tel_id] = 0.;
        }
        else
        {
            if ( fabs(syy) < 1e-8*fabs(sxx) )
                syy = 0.;
            else if ( fabs(sxx) < 1e-8*fabs(syy) )
                sxx = 0.;
            if ( sxx > syy && syy >= 0. )
            {
                image->length[tel_id] = sqrt(sxx);
                image->width[tel_id] = sqrt(syy);
                b = 0.;
                a = sy;
            }
            else if ( syy >= 0. && sxx >= 0. )
            {
                image->length[tel_id] = sqrt(syy);
                image->width[tel_id] = sqrt(sxx);
                    b = 100000.;
                    a = sy - b*sx;
            }
        }
        image->SetTelAlpha(tel_id, atan(b));  
        image->SetTelSize(tel_id, sA);
        image->SetTelImageX(tel_id, sx) ;
        image->SetTelImageY(tel_id, sy) ;
    }


}