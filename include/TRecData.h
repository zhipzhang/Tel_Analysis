#include "TTree.h"
#include "TImage_Parameter.h"
#include "TCuts.h"
#include "TMcData.h"
#include "LACTree.h"

class TRecData :public TMcData
{
public:

    TTree* rec_tree;
    float true_direction[2]; //[0] azimuth [1] altitude
    float rec_direction[2];
    bool  exist_lookup;
    float rec_energy;
    float rec_core[2];
     
    void InitWrite();
    void InitRead(TTree* t);
    void RecShower(TImage_Parameter* image, TCuts* cut_option);
    double compute_direction_error();
    TRecData(LACTree*);
    ~TRecData();



};
        
        