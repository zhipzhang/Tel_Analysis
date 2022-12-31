#include "TMcData.h"

TMcData::TMcData()
{
    primary_id = -1;
    energy = 0.;
    std::fill(core_pos,core_pos + 2, 0.);
    runnumber = 0;
    eventnumber = 0;
    memset(Tel_direction, 0, 2*LACT_MAXTEL*sizeof(double));
    memset(Tel_position, 0, 3*LACT_MAXTEL*sizeof(double));
}
void TMcData::GetData(LACTree* DSTTree)
{
    SetPrimaryID(DSTTree->primary);
    SetEnergy(DSTTree->energy);
    SetCorePos(DSTTree->xcore, DSTTree->ycore);
    SetTrueDirection(DSTTree->az, 90 - DSTTree->ze);
    SetRunNumber(DSTTree->runNumber);
    SetEventNumber(DSTTree->eventNumber);
    SetPointDirection(DSTTree->az,  90 - DSTTree->ze);
    SetWeight(DSTTree->weight);
    SetHeight(DSTTree->hmax, DSTTree->xmax, DSTTree->emax, DSTTree->cmax);
}

void TMcData::Reset()
{
    primary_id = -1;
    energy = 0.;
    std::fill(core_pos,core_pos + 2, 0.);
    runnumber = 0;
    eventnumber = 0;
    hmax = xmax = emax = cmax = 0.;
    memset(Tel_direction, 0, 2*LACT_MAXTEL*sizeof(Tel_direction[0][0]));
}
TMcData::~TMcData()
{

}