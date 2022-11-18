/*
    define the method in LACTree.h
*/

#include "LACTree.h"
LACTree::LACTree()
{
    Mctree = 0;
    EventTree = 0;
    LTrig = 0;
    Ntrig = 0;
    ntel = 0;
    ntel_data = 0;
    runNumber = 0;
    eventNumber = 0;
    primary = 0;
    energy = 0;
    xcore = 0;
    ycore = 0;
    flag = 0;
    fadc_read_write = 0;
    fillPeLeaf = 0;
    resetDataVectors();
}
bool LACTree::initMctree()
{
    Mctree = new TTree( "mc", "mc events" );
    Mctree->SetMaxTreeSize( 1000 * Long64_t( 2000000000 ) );
    Mctree->Branch( "runNumber", &runNumber, "runNumber/i" );
    Mctree->Branch( "eventNumber", &eventNumber, "eventNumber/i" );
    Mctree->Branch( "MCprim", &primary, "MCprimary/s" );
    Mctree->Branch( "MCe0", &energy, "MCenergy/F" );
    Mctree->Branch( "MCxcore", &xcore, "MCxcore/F" );
    Mctree->Branch( "MCycore", &ycore, "MCycore/F" );
    Mctree->Branch( "MCze", &ze, "MCze/F" );
    Mctree->Branch( "MCaz", &az, "MCaz/F" );
    
    return true;
}
bool LACTree::initEventTree()
{
    EventTree = new TTree("event","event data tree");
    EventTree->SetMaxTreeSize(1000 * Long64_t(2000000000));
    EventTree->Branch("ntel_data", &ntel_data, "ntel_data/i") ;
    EventTree->Branch("tel_data", tel_data,"tel_data[ntel_data]/i");
    EventTree->Branch("ntel ", &ntel,"ntel/i");
    EventTree->Branch("Paz", Point_Az, "Paz[ntel]/F");
    EventTree->Branch("Pal", Point_Al, "Pal[ntel]/F");
    char tname[100];
    EventTree->Branch("ltrig", &LTrig,"ltrig/i");
    EventTree->Branch("ntrig", &Ntrig,"ntrig/i");
    EventTree->Branch("ltrig_list", LTrig_list,"ltrig_list[ntrig]/i");
    EventTree->Branch( "runNumber", &runNumber, "runNumber/i" );
    EventTree->Branch( "eventNumber", &eventNumber, "eventNumber/i" );
    EventTree->Branch( "MCprim", &primary, "MCprimary/s" );
    EventTree->Branch( "MCe0", &energy, "MCenergy/F" );
    EventTree->Branch( "MCxcore", &xcore, "MCxcore/F" );
    EventTree->Branch( "MCycore", &ycore, "MCycore/F" );
    EventTree->Branch( "MCze", &ze, "MCze/F" );
    EventTree->Branch( "MCaz", &az, "MCaz/F" );
    EventTree->Branch( "weight", &weight,"weight/F");
    EventTree->Branch("width", width, "width[ntrig]/F");
    EventTree->Branch("length", length, "length[ntrig]/F");
    EventTree->Branch("imgx", x_img, "imgx[ntrig]/F");
    EventTree->Branch("imgy", y_img, "imgy[ntrig]/F");
    EventTree->Branch("size", size, "size[ntrig]/F");
    if(fillPeLeaf)
    {
        sprintf(tname,"Pe[ntel_data][%d]/s",LACT_MAXPIXELS);
        EventTree->Branch("Pe",pe_list,tname);
    }
    EventTree->Branch( "Write_fadc", &fadc_read_write, "Write_fadc/B");
    if(fadc_read_write)
    {
        sprintf(tname,"fadc_HG[ntel_data][%d]/S",LACT_MAXPIXELS);
        EventTree->Branch("fadc_HG", fadc_HG, tname);
        sprintf(tname,"fadc_sum[ntel_data][%d]/F",LACT_MAXPIXELS);
        EventTree->Branch("fadc_sum", fadc_sums, tname);
        sprintf(tname,"fadc_num_samples[ntel_data]/S");
        EventTree->Branch("fadc_num_samples", fadc_num_samples, tname);
        sprintf(tname,"fadc_trace[ntel_data][%d][%d]/S",LACT_SUMWINDOW,LACT_MAXPIXELS);
        EventTree->Branch("fadc_traces", fadc_trace, tname);
        sprintf(tname,"fadc_pedestal[ntel_data][%d]/F",LACT_MAXPIXELS);
        EventTree->Branch("fadc_pedestal", fadc_pedestal, tname);
    }

    
    resetDataVectors();
    return true;
}
void LACTree::resetDataVectors(unsigned int iNMAX_TEL, unsigned int iNMAX_PIXELS)
{   
    if(iNMAX_TEL>=LACT_MAXTEL)
    {
        iNMAX_TEL=LACT_MAXTEL;
    }
    if(iNMAX_PIXELS>=LACT_MAXPIXELS)
    {
        iNMAX_PIXELS=LACT_MAXPIXELS;
    }


    LTrig = 0;
    Ntrig = 0;
    for(unsigned int i=0; i<iNMAX_TEL; i++)
    {
        LTrig_list[i] = 0;
        LTime[i] = 0.;
        Point_Al[i] = 0.;
        Point_Az[i] = 0.;
    }
    memset(pe_list, 0, LACT_MAXTEL*LACT_MAXPIXELS*sizeof(pe_list[0][0]));
    memset(fadc_pedestal, 0., LACT_MAXPIXELS*LACT_MAXTEL*sizeof(fadc_pedestal[0][0]));
    memset(fadc_sums, 0., LACT_MAXTEL * LACT_MAXPIXELS * sizeof(fadc_sums[0][0]));
    memset(fadc_num_samples, 0, LACT_MAXTEL* sizeof(fadc_num_samples[0]));
    memset(fadc_trace, 0, LACT_MAXTEL * LACT_SUMWINDOW *LACT_MAXPIXELS * sizeof(fadc_trace[0][0][0]));
}

bool LACTree::initEventTree(TTree *t)
{
    EventTree = t;
    if(EventTree && EventTree->GetEntries() == 0)
    {
        std::cout << "Event Tree Have no Entries" << std::endl;
    }

    if(EventTree->GetBranch("MCe0"))
    {
        fMC = true;
    }

    EventTree->SetBranchAddress("runNumber", &runNumber);
    EventTree->SetBranchAddress("eventNumber", &eventNumber);
    //EventTree->SetBranchAddress("ntel", &ntel);
    
    EventTree->SetBranchAddress("Paz", Point_Az);
    EventTree->SetBranchAddress("Pal", Point_Al);
    EventTree->SetBranchAddress("ltrig", &LTrig);
    EventTree->SetBranchAddress("ntrig", &Ntrig);
    EventTree->SetBranchAddress("ltrig_list", &LTrig_list);
    EventTree->SetBranchAddress("ntel_data", &ntel_data);
    EventTree->SetBranchAddress("tel_data", tel_data);
    EventTree->SetBranchAddress("Write_fadc", &fadc_read_write);
    if(EventTree->GetBranch("Pe"))
    {
        EventTree->SetBranchAddress("Pe", pe_list);
        setFillPEleaf(true);
    }
    
    if( fMC )
    {
        EventTree->SetBranchAddress( "MCprim", &primary );
        EventTree->SetBranchAddress( "MCe0", &energy );
        EventTree->SetBranchAddress( "MCxcore", &xcore );
        EventTree->SetBranchAddress( "MCycore", &ycore );
        EventTree->SetBranchAddress( "MCze", &ze );
        EventTree->SetBranchAddress( "MCaz", &az );
    }



}
