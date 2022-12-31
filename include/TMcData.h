#ifndef _MC_DATA_
#define _MC_DATA_
#include "LACTree.h"
#include "Limits_defined.h"
#include "TMath.h"
class TMcData
{
    public:
        int primary_id;
        int ntel;
        double energy;
        double core_pos[2];
        double true_direction[2];
        double point_direction[2];
        double Tel_direction[LACT_MAXTEL][2];
        double Tel_position[LACT_MAXTEL][3];
        double Rp[LACT_MAXTEL];
        double weight;
        int runnumber;
        int eventnumber;
        double hmax;
        double xmax;
        double cmax;
        double emax;

        TMcData();
        ~TMcData();
        void GetData(LACTree* DSTTree);
        void SetPrimaryID(int id)
        {
            primary_id = id;
        }
        void SetEnergy(double ienergy)
        {
            energy = ienergy;
        }
        void SetCorePos(double corex, double corey)
        {
            core_pos[0] = corex;
            core_pos[1] = corey;
        }
        void SetTrueDirection(double azimuth, double altitude)
        {
            true_direction[0] = azimuth;
            true_direction[1] = altitude;
        }
        void SetRunNumber(int run)
        {
            runnumber = run;
        }
        void SetEventNumber(int ievent)
        {
            eventnumber = ievent;
        }
        void SetPointDirection(double point_az, double point_al)
        {
            point_direction[0] = point_az ;
            point_direction[1] = point_al ;
        }
        void SetTelPos(double tel_pos[][3])
        {
            memcpy(Tel_position, tel_pos, 3 * LACT_MAXTEL* sizeof(tel_pos[0][0]));
        }
        void SetTelDirection(double tel_dir[][2])
        {
            memcpy(Tel_direction, tel_dir, 2 * LACT_MAXTEL * sizeof(tel_dir[0][0]));
        }
        void SetWeight(double w)
        {
            weight = w;
        }
        void SetHeight(double shower_hmax, double shower_xmax, double shower_emax, double shower_cmax)
        {
            hmax = shower_hmax;
            xmax = shower_xmax;
            emax = shower_emax;
            cmax = shower_cmax;
        }
        
        int GetRunNumber()
        {
            return runnumber;
        }
        int GetEventNumber()
        {
            return eventnumber;
        }
        double GetEnergy()
        {
            return energy;
        }
        double GetAzimuth()
        {
            return true_direction[0];
        }
        double GetAltitude()
        {
            return true_direction[1];
        }
        void Reset();
        


};
#endif