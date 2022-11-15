#ifndef _REC_DATA_
#define _REC_DATA_

#include "TTree.h"
#include "TImage_Parameter.h"
#include "TCuts.h"
#include "TMcData.h"
#include "LACTree.h"

class TRecData :public TMcData
{
public:

    TTree* rec_tree;
    //double true_direction[2]; //[0] azimuth [1] altitude
    double rec_direction[2];
    bool  exist_lookup;
    double rec_energy;
    double rec_core[2];
     
    void InitWrite();
    void InitRead(TTree* t);
    void RecShower(TImage_Parameter* image, TCuts* cut_option);
    TTree* GetRecTree()
    {
        return rec_tree;
    }
    double compute_direction_error();
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
    void Reset();
    TRecData();
    ~TRecData();



};
        
        
#endif