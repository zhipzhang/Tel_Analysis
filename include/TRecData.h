#ifndef _REC_DATA_
#define _REC_DATA_

#include "TTree.h"
#include "TImage_Parameter.h"
#include "TUserCuts.h"
#include "TMcData.h"
#include "LACTree.h"

class TRecData :public TMcData
{
public:

    TTree* rec_tree;
    //double true_direction[2]; //[0] azimuth [1] altitude
    double rec_direction[2];
    int   exist_lookup;
    double rec_energy;
    double rec_core[2];
    double direction_error;
    double theta2;
    double sourcex;
    double sourcey;
    double beta;
    double rec_x; //intersection points coordinates
    double rec_y;

    double MRSL;
    double MRSW;
     
    void InitWrite();
    void InitRead(TTree* t);
    bool RecShower(TImage_Parameter* image, TUserCuts* cut_option);
    TTree* GetRecTree()
    {
        return rec_tree;
    }
    void compute_direction_error();
    void GetData(LACTree*);
    void SetRecDirection(double ref_az, double ref_al)
    {
        rec_direction[0] = ref_az;
        rec_direction[1] = ref_al;
    }
    void SetRecCore(double x, double y)
    {
        rec_core[0] = x;
        rec_core[1] = y;
    }
    void SetLookup()
    {
        exist_lookup = 1;
    }
    void SetRecEnergy( double rec)
    {
        rec_energy = rec;
    }
    void SetShape(double w, double l)
    {
        MRSL = l, MRSW = w;
    }
    void SetRecPoint(double x, double y)
    {
        rec_x = x, rec_y = y;
    }
    void ReWeight( double *, double *);
    
    double GetDirectionError()
    {
        return direction_error;
    }
    double GetCorex()
    {
        return core_pos[0];
    }
    double GetCorey()
    {
        return core_pos[1];
    }
    double GetCoredist()
    {
        return sqrt(pow(core_pos[0], 2) + pow(core_pos[1], 2));
    }
    void Reset();
    TRecData();
    ~TRecData();



};
        
        
#endif