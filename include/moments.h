#ifndef _moments_
#define _moments_
#include "TImage_Parameter.h"
#include "Limits_defined.h"
#include <map>
#include <vector>
void image_clean(double pe_list[][LACT_MAXPIXELS], std::map<int, std::vector<int> > *pixel_neighbors, double* tail_cuts, 
                                                std::vector<int> *pixel_in_image, unsigned int ntel, unsigned int* Trig_List, int* npix);
void compute_moments(TImage_Parameter* image,double pe_list[][LACT_MAXPIXELS], std::vector<int>* pixel_in_image, unsigned int ntel, 
                        unsigned int* Trig_List, double x_pix[][LACT_MAXPIXELS], double y_pix[][LACT_MAXPIXELS]);






#endif