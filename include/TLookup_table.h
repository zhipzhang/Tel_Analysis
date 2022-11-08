#include "TTree.h"
#include "TH2D.h"
#include "TImage_Parameter.h"
#include <vector>
#include "TRecData.h"



class TLookup_table
{
    public:
        TTree* lookup_tree;
        std::vector<TH2D*> energy_lookup;
        std::vector<TH2D*> width_lookup;
        std::vector<TH2D*> length_lookup;
        std::vector<double> zenith_distribution;
        

        void Init();
        void Init(TTree* t);
        void GenerateLookup(TTree* image_tree, TRecData* rec);


        TLookup_table();
        ~TLookup_table();

};