#ifndef _moments_
#define _moments_
#include "TImage_Parameter.h"
#include "Limits_defined.h"
#include "TRecData.h"
#include <map>
#include <vector>
#include <queue>
void image_clean(double pe_list[][LACT_MAXPIXELS], std::map<int, std::vector<int> > *pixel_neighbors, double* tail_cuts, 
                                                std::vector<int> *pixel_in_image, unsigned int ntel, unsigned int* Trig_List, int* npix);
void image_clean2(double pe_list[][LACT_MAXPIXELS], std::map<int, std::vector<int> > *pixel_neighbors, double* tail_cuts, 
                                                std::vector<int> *pixel_in_image, unsigned int ntel, unsigned int* Trig_List, int* npix);
void compute_moments(TImage_Parameter* image, TRecData* rec, double pe_list[][LACT_MAXPIXELS], std::vector<int>* pixel_in_image, unsigned int ntel, 
                        unsigned int* Trig_List, float x_pix[][LACT_MAXPIXELS], float y_pix[][LACT_MAXPIXELS]);
void ComputeRp(TImage_Parameter* image, TRecData* rec, double tel_pos[][3]);





#endif