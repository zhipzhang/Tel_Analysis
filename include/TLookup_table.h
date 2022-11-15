#include "TTree.h"
#include "TH2D.h"
#include "TImage_Parameter.h"
#include <vector>
#include "TRecData.h"



class TLookup_table
{
    public:
        TTree* lookup_tree;
        TH2D energy_lookup;
        TH2D width_lookup;
        TH2D length_lookup;
        std::vector<double> zenith_distribution;
        

        void InitHist();
        void GenerateLookup(TTree* image_tree, TRecData* rec);


        TLookup_table();
        ~TLookup_table();

};