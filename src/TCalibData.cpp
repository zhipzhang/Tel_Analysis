#include "TCalibData.h"
#include <iostream>


TCalibData::TCalibData()
{
    Tel_id = -1;
    Npix = 0;
    Num_sumwindow = 0;
    Low_gain = 1;
    memset(High_pedstal, 0, LACT_MAXPIXELS*sizeof(High_pedstal[0]));
    memset(High_pedstal_var, 0, LACT_MAXPIXELS*sizeof(High_pedstal_var[0]));
    memset(Low_pedstal, 0, LACT_MAXPIXELS*sizeof(Low_pedstal[0]));
    memset(Low_pedstal_var, 0, LACT_MAXPIXELS*sizeof(Low_pedstal_var[0]));
    memset(Low_Convert, 0, LACT_MAXPIXELS*sizeof(Low_Convert[0]));
    memset(High_Convert, 0, LACT_MAXPIXELS*sizeof(High_Convert[0]));

}

void TCalibData::InitWrite()
{
    calibration_tree = new TTree("calibration_data", "calib_data");
    calibration_tree->Branch("Tel_id", &Tel_id, "Tel_id/I");
    calibration_tree->Branch("Npix", &Npix, "Npix/I");
   // calibration_tree->Branch("num_sumwindow", &Num_sumwindow, "num_sumwindow/I");
    calibration_tree->Branch("High_pedstal", High_pedstal, "High_pedestal[Npix]/F");
    calibration_tree->Branch("High_pedstal_var", High_pedstal_var, "High_pedestal_var[Npix]/F");
    calibration_tree->Branch("High_Convert", High_Convert, "High_Convert[Npix]/F");
    calibration_tree->Branch("Low_Convert", Low_Convert, "Low_Convert[Npix]/F");
    calibration_tree->Branch("Low_pedestal", Low_pedstal, "Low_pedestal[Npix]/F");
    calibration_tree->Branch("Low_pedestal_var", Low_pedstal_var, "Low_pedestal_var[Npix]/F");



}

void TCalibData::FillData(AllHessData *hsdata)
{
    for(int itel = 0; itel < hsdata->run_header.ntel; itel++)
    {
        Tel_id = hsdata->tel_moni[itel].tel_id;
        Npix = hsdata->tel_moni[itel].num_pixels;
        if(Npix > LACT_MAXPIXELS)
        {
            std::cerr << "The number of Pixels is too Large ! " << std::endl;
            exit(EXIT_FAILURE);
        }
        for(unsigned int p = 0; p < Npix; p++)
        {
            if(hsdata->tel_moni[itel].num_ped_slices > 0)
            {
                High_pedstal[p] = hsdata->tel_moni[itel].pedestal[HI_GAIN][p] / (double) hsdata->tel_moni[itel].num_ped_slices;
                Low_pedstal[p] = hsdata->tel_moni[itel].pedestal[LOW_GAIN][p] /(double) hsdata->tel_moni[itel].num_ped_slices;
            }
            else
            {
                High_pedstal[p] = 0;
                Low_pedstal[p] = 0;
            }
            High_pedstal_var[p] = hsdata->tel_moni[itel].noise[HI_GAIN][p];
            Low_pedstal_var[p]  = hsdata->tel_moni[itel].noise[LOW_GAIN][p];
            High_Convert[p] = hsdata->tel_lascal[itel].calib[HI_GAIN][p] * CALIB_SCALE;
            Low_Convert[p] = hsdata->tel_lascal[itel].calib[LOW_GAIN][p] * CALIB_SCALE;
        }
        calibration_tree->Fill();

    }
}