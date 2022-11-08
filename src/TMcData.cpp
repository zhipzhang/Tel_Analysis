#include "TMcData.h"

TMcData::TMcData(LACTree* DSTTree)
{
    SetPrimaryID(DSTTree->primary);
    SetEnergy(DSTTree->energy);
    SetCorePos(DSTTree->xcore, DSTTree->ycore);
    SetTrueDirection(DSTTree->az, DSTTree->altitude);
    SetRunNumber(DSTTree->runNumber);
    SetEventNumber(DSTTree->eventNumber);
    SetPointDirection(DSTTree->Tel_ze, DSTTree->Tel_al);
}