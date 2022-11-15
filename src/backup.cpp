bool DST_fillEvent(LACTree* fData, AllHessData* hsdata);

bool DST_fillEvent(LACTree* fData, AllHessData* hsdata)
{
    if(!fData || !hsdata)
    {
        return false;
    }

    fData->ntel = (unsigned short int)hsdata->run_header.ntel;
    fData->resetDataVectors(LACT_MAXTEL,LACT_MAXPIXELS);
    fData->eventNumber = hsdata->mc_event.event;
    fData->runNumber = hsdata->run_header.run;
 
    //only write MC data frist 
    fData->primary = hsdata->mc_shower.primary_id;
    fData->energy = hsdata->mc_shower.energy;
    fData->az = hsdata->mc_shower.azimuth * TMath::RadToDeg();
    fData->ze = 90. - hsdata->mc_shower.altitude * TMath::RadToDeg(); 
    fData->xcore = hsdata->mc_event.xcore;
    fData->ycore = hsdata->mc_event.ycore;  
    double lg_E = log10(fData->energy);
    double index = hsdata->mc_run_header.spectral_index;
    fData->weight = pow(fData->energy, (-2.7 - index));

    //check flag whether it is trigger;
    //trigger data

    fData->Ntrig = hsdata->event.central.num_teltrg;
    unsigned int i_ntel_trig = 0;
    bitset<8 * sizeof(unsigned long)> i_localTrigger;//use bit to represent the trigger_list?

    //loop over all trigger event(not write now)
    for(unsigned int t=0; t<(unsigned int )hsdata->event.central.num_teltrg; t++)
    {
        if(hsdata->event.central.teltrg_list[t] < (int)i_localTrigger.size())   
        {
            i_localTrigger.set(hsdata->event.central.teltrg_list[t] - 1, true);
        }
        if(t < (unsigned int)hsdata->run_header.ntel)
        {
            fData->LTrig_list[i_ntel_trig] = hsdata->event.central.teltrg_list[t];
            fData->LTime[i_ntel_trig] = hsdata->event.central.teltrg_time[t];
        }
        i_ntel_trig++;
    }

    // tracking data
    unsigned int z =0;
    for( unsigned int i=0; i<(unsigned short int)hsdata->run_header.ntel;i++)
    {
        fData->Point_Az[z] = hsdata->event.trackdata[i].azimuth_raw * TMath::RadToDeg();
        fData->Point_Al[z] = hsdata->event.trackdata[i].altitude_raw * TMath::RadToDeg();
        z++;
    }

    //event data 
    fData->ntel_data = (unsigned short int)hsdata->event.central.num_teldata;
    unsigned short int i_ntel_data = 0;

    assert(fData->ntel_data > 0 );
    //loop over all telescopes
    for( unsigned short int i = 0; i<fData->ntel_data; i++)
    {
        fData->tel_data[i_ntel_data] = (unsigned short int)hsdata->event.central.teldata_list[i];
        unsigned int telID = fData->tel_data[i_ntel_data]-1;
        
        fData->fadc_num_samples[i_ntel_data]  = (unsigned short int)hsdata->event.teldata[telID].raw->num_samples;
        fData->Telescope_ZeroSuppression[i_ntel_data] = (unsigned short int)hsdata->event.teldata[telID].raw->zero_sup_mode;

        if( hsdata->camera_set[telID].num_pixels >= LACT_MAXPIXELS)
        {
            cout << "Error: There are too much pixels in Tel" << telID << endl;
            cout << "The limit is "<< LACT_MAXTEL << " but there are " << hsdata->camera_set[telID].num_pixels << endl;  
        }
        // loop all pixels
        for(int p=0; p<hsdata->camera_set[telID].num_pixels; p++)
        {
            if(fFillPELeaf)
            {
                fData->pe_list[i_ntel_data][p] = hsdata->mc_event.mc_pe_list[telID].pe_count[p];
            }
            if(!NO_fadc)
            {
                fData->fadc_HG[i_ntel_data][p] = (unsigned short int) (hsdata->event.teldata[telID].raw->adc_sum[HI_GAIN][p]);

                if( hsdata->event.teldata[telID].raw->adc_sample && hsdata->event.teldata[telID].raw->adc_sample[HI_GAIN] 
                        && hsdata->event.teldata[telID].raw->num_samples)
                {
                    if( hsdata->event.teldata[telID].raw->num_gains == 2)
                    {
                        Only_High_Gain = false;
                    }

                    // Now only treat the High Gain Situation , May Add Low_Gain next time
                        for( int  t = 0; t < hsdata->event.teldata[telID].raw->num_samples; t++)
                        {
                            fData->fadc_trace[i_ntel_data][t][p] = hsdata->event.teldata[telID].raw->adc_sample[HI_GAIN][p][t];
                        }

                }
                else // no trace data
                {
                    fData->fadc_sums[i_ntel_data][p] = hsdata->event.teldata[telID].raw->adc_sum[HI_GAIN][p] - hsdata->tel_moni[telID].pedestal[HI_GAIN][p];
                }

                //record pedestal 
                if(hsdata->tel_moni[telID].num_ped_slices > 0)
                {
                    fData->fadc_pedestal[i_ntel_data][p] = hsdata->tel_moni->pedestal[HI_GAIN][p]/(double) (hsdata->tel_moni[telID].num_ped_slices);
                }
            }
        }
        
        /*
        fData->Rp[i_ntel_data] = line_point_distance(fData->xcore, fData->ycore, 0.,cos(hsdata->mc_shower.altitude)*
        cos(hsdata->mc_shower.azimuth), cos(hsdata->mc_shower.altitude) * sin(-1 * hsdata->mc_shower.azimuth), sin(hsdata->mc_shower.altitude),
        hsdata->run_header.tel_pos[telID][0], hsdata->run_header.tel_pos[telID][1], hsdata->run_header.tel_pos[telID][2]);
        double theta = compute_true_dsp(hsdata->mc_shower.azimuth, hsdata->mc_shower.altitude, 
        fData->xcore, fData->ycore,hsdata->run_header.tel_pos[telID][0], hsdata->run_header.tel_pos[telID][1],
        hsdata->event.trackdata[telID].azimuth_raw, hsdata->event.trackdata[telID].altitude_raw);
        fData->dsp_error[i_ntel_data] = compute_axis_diff(theta, hsdata->event.teldata[telID].img[0].phi);
        */
        i_ntel_data++;
    }    
