#include "moments.h"
#include <algorithm>
#include <cmath>
#include <queue>
#include "TMath.h"
#include "TMathBase.h"
#include "rec_tools.h"
#include "TRecData.h"

void image_clean(double pe_list[][LACT_MAXPIXELS], std::map<int, std::vector<int> > *pixel_neighbors, double* tail_cuts, 
                                                std::vector<int> *pixel_in_image, unsigned int ntel, unsigned int* Trig_List, int* npix)
{
    unsigned int tel_id;
    for(int i = 0; i < ntel; i++)
    {
        tel_id = Trig_List[i] - 1;
        for(int p = 0; p < npix[tel_id]; p++)
        {
            if( pe_list[tel_id][p] > tail_cuts[1])
            {
                for(auto k : pixel_neighbors[tel_id][p])
                {
                    if(pe_list[tel_id][k] > tail_cuts[0])
                    {
                        pixel_in_image[tel_id].push_back(p);
                        break;
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
                        break;
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


void image_clean2(double pe_list[][LACT_MAXPIXELS], std::map<int, std::vector<int> > *pixel_neighbors, double* tail_cuts, 
                                                std::vector<int> *pixel_in_image, unsigned int ntel, unsigned int* Trig_List, int* npix)
{
    unsigned int tel_id;
    bool flag = 1;
    for(int i = 0; i < ntel; i++)
    {
        tel_id = Trig_List[i] - 1;
        std::queue<int> tmp;
        bool *res = (bool*)calloc(npix[tel_id], sizeof(bool));
        int max_pos = std::max_element(pe_list[tel_id], pe_list[tel_id] + npix[tel_id]) - pe_list[tel_id];
        if( pe_list[tel_id][max_pos] < 10)
        {
            continue;
        }
        else {
            tmp.push(max_pos);
        }
        while( !tmp.empty() )
        {
            int tmp_p = tmp.front();
            res[tmp_p] = 1;
            pixel_in_image[tel_id].push_back(tmp_p);
            tmp.pop();
            if( pe_list[tel_id][tmp_p] > tail_cuts[1])
            {
                for(auto k : pixel_neighbors[tel_id][tmp_p])
                {
                    if(res[k])
                    {
                        continue;
                    }
                    else
                    {
                        if(pe_list[tel_id][k] > tail_cuts[0])
                        {
                            res[k] = 1;
                            tmp.push(k);
                        }
                    }

               }
            }
            else if(pe_list[tel_id][tmp_p] > tail_cuts[0])
           {
                for(auto k : pixel_neighbors[tel_id][tmp_p])
                {
                    if(res[k])
                    {
                        continue;
                    }
                    else
                    {
                        if(pe_list[tel_id][k] > tail_cuts[1])
                        {
                            res[k] = 1;
                            tmp.push(k);
                        }
                    }

               }

           }
       }
       free(res);
    }

}


/*be careful alpha not consider the camera rotation !*/

void compute_moments(TImage_Parameter* image, TRecData* rec, double pe_list[][LACT_MAXPIXELS], std::vector<int>* pixel_in_image,unsigned int ntel, unsigned int* Trig_List,
                     float x_pix[][LACT_MAXPIXELS], float y_pix[][LACT_MAXPIXELS], bool* known , int num_only)
{
    int tel_id;
    for(int i = 0; i < ntel; i++)
    {
        double sx = 0, sxx = 0, sxy = 0, sy = 0, syy  = 0, sA = 0. ,sx3 = 0.;
        double beta, cb, sb;
        double a = 0, b = 0;
        tel_id = Trig_List[i] - 1;
        if( num_only >0 )
        {
            if(!known[tel_id])
            {
                continue;
            }
        }
        if(pixel_in_image[tel_id].size() < 3 )
        {
            continue;
        }
        else
        {
            image->image_tel.push_back(tel_id);
        }
        double dx,dy;
        angles_to_offset(rec->GetAzimuth() * TMath::DegToRad(), rec->GetAltitude() * TMath::DegToRad(), rec->Tel_direction[tel_id][0] * TMath::DegToRad()
        , rec->Tel_direction[tel_id][1] * TMath::DegToRad(), image->GetTelFocal(tel_id) , &dx, &dy);
        image->SetTelSource(tel_id, dx, dy);
        for(int p = 0; p < pixel_in_image[tel_id].size(); p++)
        {
            int ipix = pixel_in_image[tel_id][p];
            double x = x_pix[tel_id][ipix] - dx;
            double y = y_pix[tel_id][ipix] - dy;
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
        if (fabs(sxy) > 1e-3 * fabs(sxx) && fabs(sxy) > 1e-3 * fabs(syy) )
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
        beta = atan(b);
        cb = cos(beta);
        sb = sin(beta);
        image->SetTelSize(tel_id, sA);
        image->SetTelImageX(tel_id, sx + dx) ;
        image->SetTelImageY(tel_id, sy + dy) ;
        sxx = sx3 = 0.;
        for(int j = 0; j < pixel_in_image[tel_id].size(); j++)
        {
            int ipix = pixel_in_image[tel_id][j];
            double x = x_pix[tel_id][ipix];
            double y = y_pix[tel_id][ipix];
            double A = pe_list[tel_id][ipix];
            double xp;
            xp = cb*(x -sx) + sb * (y - sy);
            sxx += (A*xp) * xp;
            sx3 += ((A*xp) * xp) * xp;
            
        }
        if(sx3 / pow(sxx, 1.5) <0. )
        {
            beta = beta * TMath::RadToDeg() + 180;
        }
        else
        {
            beta = beta * TMath::RadToDeg();
        }
        image->SetTelAlpha(tel_id, beta*TMath::DegToRad());  
    }
    
    
}

void ComputeRp(TImage_Parameter* image, TRecData* rec, double tel_pos[][3])
{

    for(int i = 0; i < image->image_tel.size(); i++)
    {
        int itel = image->image_tel[i];
        image->rp[itel] = line_point_distance(rec->core_pos[0], rec->core_pos[1], 0, cos(rec->true_direction[0] * TMath::DegToRad()) * cos(rec->true_direction[1] * TMath::DegToRad()),
                                        -sin(rec->true_direction[0] * TMath::DegToRad()) * cos(rec->true_direction[1] * TMath::DegToRad()), sin(rec->true_direction[1] * TMath::DegToRad() )
                                        ,tel_pos[itel][0], tel_pos[itel][1], tel_pos[itel][2]);
         
        image->rec_rp[itel] = line_point_distance(rec->rec_core[0], rec->rec_core[1], 0, cos(rec->rec_direction[0] * TMath::DegToRad()) * cos(rec->rec_direction[1] * TMath::DegToRad()), 
                                        -sin(rec->rec_direction[0] * TMath::DegToRad()) * cos(rec->rec_direction[1] * TMath::DegToRad()), sin(rec->rec_direction[1] * TMath::DegToRad())
                                        , tel_pos[itel][0], tel_pos[itel][1], tel_pos[itel][2]);
    }
}