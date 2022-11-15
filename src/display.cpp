#include "TH2Poly.h"
#include "TImage_Parameter.h"
#include "vector"
#include "TTree.h"
#include "TCanvas.h"
#include "TRecData.h"
#include "Limits_defined.h"
#include "TEllipse.h"
#include "TPaveText.h"

void display(TImage_Parameter* , double pe[][LACT_MAXPIXELS], std::vector<int> *pixel_in_image, double xpix[][LACT_MAXPIXELS], double ypix[][LACT_MAXPIXELS],double pix_size[][LACT_MAXPIXELS], TRecData* rec, int id);

void display(TImage_Parameter* image, double pe[][LACT_MAXPIXELS], std::vector<int> *pixel_in_image, double xpix[][LACT_MAX_TIMELEVELS], double ypix[][LACT_MAXPIXELS], double pix_size[][LACT_MAXPIXELS], TRecData* rec, int id)
{
    int num = 0;
    for(int i = 0; i < image->image_tel.size(); i++)
    {
        int tel_id = image->image_tel[i];
        if(image->GetTelSize(tel_id) > 50)
        {
            TCanvas* camera_image = new TCanvas(Form("LACT image of camera %d", tel_id), "LACT image", 1600, 1600);
            TH2Poly* camera = new TH2Poly(Form("camera %d", tel_id), "", -6, 6, -6, 6);
            for(auto j: pixel_in_image[tel_id])
            {
                double binsize = pix_size[tel_id][j] / image->GetTelFocal(tel_id);
                double x = xpix[tel_id][j];
                double y = ypix[tel_id][j];
                double bin_x[4] = {x - 0.5 * binsize, x + 0.5*binsize, x + 0.5*binsize, x - 0.5*binsize};
                double bin_y[4] = {y - 0.5 * binsize, y - 0.5*binsize, y + 0.5*binsize, y + 0.5*binsize};
                camera->AddBin(4, bin_x, bin_y);
                camera->Fill(x, y, pe[tel_id][j]);

            }
            camera->Draw("CLOZ");
            TEllipse* ellipse = new TEllipse(image->GetTelImageX(tel_id), image->GetTelImageY(tel_id), image->GetTelLength(tel_id), image->GetTelwidth(tel_id),
                                            0, 360, image->GetTelAlpha(tel_id));
            ellipse->SetLineWidth(2);
            ellipse->SetLineColor(2);
            ellipse->SetFillStyle(0);
            ellipse->Draw();

            TPaveText *pavet = new TPaveText(-6, 6.3, 6, 7.6);
            pavet->SetFillStyle(0);
            pavet->AddText(Form("event_number: %d, Tel: %d, energy: %.4lf , azimuth: %.4lf, zenith: %.4lf Rp:%.4lf m", rec->GetEventNumber(), rec->GetEnergy(), rec->GetAzimuth(), rec->GetAltitude(), image->GetTelRp(tel_id)));
            pavet->Draw("same");
            camera_image->SaveAs(Form("./image_camera%d.png", tel_id));
        }

    }
}