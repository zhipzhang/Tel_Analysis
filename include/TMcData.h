#ifndef _MC_DATA_
#define _MC_DATA_
#include "LACTree.h"
#include "Limits_defined.h"
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
        int runnumber;
        int eventnumber;

        TMcData(LACTree* DSTTree);
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
            point_direction[0] = point_az;
            point_direction[1] = point_al;
        }
        


};
#endif