#ifndef _HEELLO_
#define _HEELLO_
#include "TH2Poly.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TPaveText.h"

void display(TImage_Parameter* img, double pe[][LACT_MAXPIXELS], std::vector<int> *pixel_in_image, float xpix[][LACT_MAXPIXELS], float ypix[][LACT_MAXPIXELS],double *pix_size, TRecData* rec );

#endif