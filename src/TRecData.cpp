#include "TRecData.h"
#include "rec_tools.h"

TRecData::TRecData(LACTree* DSTTree) : TMcData(DSTTree)
{
    std::fill(rec_core, rec_core + 2, 0.);
    std::fill(rec_direction, rec_direction + 2, 0.);
    rec_energy = 0.;

}
void TRecData::InitWrite()
{
    rec_tree->Branch("true_direction", true_direction);
    rec_tree->Branch("rec_direction", rec_direction);
    rec_tree->Branch("true_energy", &energy);
    if(exist_lookup)
        rec_tree->Branch("rec_energy", &rec_energy);
    rec_tree->Branch("true_core", &core_pos);
    rec_tree->Branch("rec_core", &rec_core);

}

void TRecData::InitRead(TTree *t)
{
    rec_tree = t;
    rec_tree->SetBranchAddress("true_direction", true_direction);
    rec_tree->SetBranchAddress("rec_direction", rec_direction);
    rec_tree->SetBranchAddress("true_energy", &energy);
    if(exist_lookup)
    {
        rec_tree->SetBranchAddress("rec_energy", &rec_energy);
    }
    rec_tree->SetBranchAddress("true_core", &core_pos);
    rec_tree->SetBranchAddress("rec_core", &rec_core);
}

double TRecData::compute_direction_error()
{
    return angle_between(true_direction[1], true_direction[0], rec_direction[1], rec_direction[0]);

}

// ! Not consider the camera rotation Now  (cam_to_ref function)


void TRecData::RecShower(TImage_Parameter* image, TCuts* cut_option)
{
    double x_ref[LACT_MAXTEL], y_ref[LACT_MAXTEL], alpha_ref[LACT_MAXTEL];
    double w ; //weight we will use
    double sum_xs, sum_xs2, sum_ys, sum_ys2, sum_w;
    double xs, ys, angs; // intersect point and the angle between two lines
    double amp_red;
    double rec_az, rec_alt;
    double trans[3][3];
    double xh, yh, zh;
    double xc, yc;
    double xt[LACT_MAXTEL], yt[LACT_MAXTEL]; // tel_pos after trans
    for(int i = 0; i < image->image_tel.size(); i++)
    {
        int itel = image->image_tel[i];
        if(image->GetTelSize(itel) < cut_option->min_size)
        {
            continue;
        }

        cam_to_ref(image->GetTelImageX(itel), image->GetTelImageY(itel), image->GetTelAlpha(itel), point_direction[0], point_direction[1],
            0., Tel_direction[itel][0], Tel_direction[itel][1], image->GetTelFocal(itel), &x_ref[itel], &y_ref[itel], &alpha_ref[itel]);
    }
    sum_xs = sum_ys = sum_w = sum_xs2 = sum_ys2 = 0.;
    for(int i = 0; i < image->image_tel.size(); i++)
    {
        int itel = image->image_tel[i];
        if(image->GetTelSize(itel) < cut_option->min_size || image->GetTelRp(itel) > cut_option->max_rp || 
                image->GetDist(itel) > cut_option->max_dist)
        {
            continue;
        }
        for(int j = 0; j < i; j++)
        {
            int jtel = image->image_tel[j];
            if( intersect_lines(x_ref[itel], y_ref[itel], alpha_ref[itel], x_ref[jtel], y_ref[jtel], alpha_ref[jtel],
                                    &xs, &ys, &angs) != 1)
            {
                continue;
            }
            amp_red = (image->GetTelSize(itel) * image->GetTelSize(jtel)) / (image->GetTelSize(itel) + image->GetTelSize(jtel));
            w = pow(amp_red * sin(angs) * (1 - image->GetTelwidth(itel)/image->GetTelLength(itel)) * (1 - image->GetTelwidth(jtel)/image->GetTelLength(jtel)), 2);
            sum_w   += w;
            sum_xs  += xs * w;
            sum_xs2 += xs * xs * w;
            sum_ys  += ys * w;
            sum_ys2 += ys * ys * w;
        }
    }
    if( fabs(sum_w) < 1e-10)
    {
        std::cerr << "Event number" << eventnumber << "of Runnumber " << runnumber << " Direction Rec Failed !" << std::endl;
        return;
    }
    sum_xs /= sum_w;
    sum_ys /= sum_w;
    offset_to_angles(sum_xs, sum_ys, point_direction[0], point_direction[1], 1.0, &rec_az, &rec_alt);
    rec_az -= (2.*TMath::Pi()) * floor(rec_az / (2 * TMath::Pi()));
    SetPointDirection(rec_az, rec_alt);
    get_shower_trans_matrix(rec_az, rec_alt, trans);
    for( int i = 0; i < image->image_tel.size(); i++)
    {
        int itel = image->image_tel[i];
        xt[itel] = trans[0][0]*Tel_position[itel][0] + 
                   trans[0][1]*Tel_position[itel][1] +
                   trans[0][2]*Tel_position[itel][2];
        yt[itel] = trans[1][0]*Tel_position[itel][0] + 
                   trans[1][1]*Tel_position[itel][1] +
                   trans[1][2]*Tel_position[itel][2];

    }

    sum_xs = sum_ys = sum_w = sum_xs2 = sum_ys2 = 0.;
    for(int i = 0; i < image->image_tel.size(); i++)
    {
        int itel = image->image_tel[i];
        for(int j = 0; j < i; j++)
        {
            int jtel = image->image_tel[j];
            if(intersect_lines(xt[itel], yt[itel], image->GetTelAlpha(itel), 
                                xt[jtel], yt[jtel], image->GetTelAlpha(jtel), &xs, &ys, &angs) != 1)
            {
                continue;
            }
            amp_red = (image->GetTelSize(itel) * image->GetTelSize(jtel)) / (image->GetTelSize(itel) + image->GetTelSize(jtel));
            w = pow(amp_red * sin(angs) * (1 - image->GetTelwidth(itel)/image->GetTelLength(itel)) * (1 - image->GetTelwidth(jtel)/image->GetTelLength(jtel)), 2);
            sum_w   += w;
            sum_xs  += xs * w;
            sum_xs2 += xs * xs * w;
            sum_ys  += ys * w;
            sum_ys2 += ys * ys * w;

        }
    }
    if( sum_w == 0)
    {
        return;
    }
    xs = sum_xs / sum_w;
    ys = sum_ys / sum_w;
    xh = trans[0][0] * xs +
        trans[1][0] * ys;
    yh = trans[0][1] * xs +
        trans[1][1] * ys;
    zh = trans[0][2] * xs +
        trans[1][2] * ys;
    
    xc = xh - trans[2][0]*zh/trans[2][2];
    yc = yh - trans[2][1]*zh/trans[2][2];
    SetCorePos(xc, yc);

}