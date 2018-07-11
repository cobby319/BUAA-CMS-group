#include "Trigger.h"
namespace trigger
{
  int trigDoubleMu[4] = {8,9,10,11};
  int trigSingleMu[4] = {10,11,15,16};
  int trigDoubleE[2] = {12,13};
  int trigHighPtE[1] = {17};//Located in DoubleElectron
  int trigSingleE[2] = {11,12};
  int trigEMu[2] = {0,3}; //all DZ paths are still missing, as well as Mu12Ele23
  int trigSinglePhoton[12] = {1,21,0,20,32,31,30,29,28,27,26,25}; //ordered by decreasing Pt, see below.

  //The pt threshold corresponding to 
  int lowThresholdPt[12] = {300,300,250,250,165,120,90,75,50,36,30,22}; //ordered by Pt to avoid double counting.
  int highThresholdPt[12] = {9999,9999,300,300,250,165,120,90,75,50,36,30}; //ordered by Pt to avoid double counting.

  int passTrigger(int trig, ULong64_t TrigHltDiMu, ULong64_t TrigHltMu, ULong64_t TrigHltDiEl, ULong64_t TrigHltEl, ULong64_t TrigHltElMu, ULong64_t TrigHltPhot, std::vector<unsigned int> *TrigHltDiMu_prescale, std::vector<unsigned int> *TrigHltMu_prescale, std::vector<unsigned int> *TrigHltDiEl_prescale, std::vector<unsigned int> *TrigHltEl_prescale, std::vector<unsigned int> *TrigHltElMu_prescale, std::vector<unsigned int> *TrigHltPhot_prescale, double selectedPhotonPt){
    std::vector<std::vector<int> > trigList(MC_Photon);
    trigList[DoubleMu].insert(trigList[DoubleMu].end(),trigDoubleMu,trigDoubleMu+(sizeof(trigDoubleMu)/sizeof(trigDoubleMu[0])));
    trigList[SingleMu].insert(trigList[SingleMu].end(),trigSingleMu,trigSingleMu+(sizeof(trigSingleMu)/sizeof(trigSingleMu[0])));
    trigList[DoubleE].insert(trigList[DoubleE].end(),trigDoubleE,trigDoubleE+(sizeof(trigDoubleE)/sizeof(trigDoubleE[0])));
    trigList[HighPtE].insert(trigList[HighPtE].end(),trigHighPtE,trigHighPtE+(sizeof(trigHighPtE)/sizeof(trigHighPtE[0])));
    trigList[SingleE].insert(trigList[SingleE].end(),trigSingleE,trigSingleE+(sizeof(trigSingleE)/sizeof(trigSingleE[0])));
    trigList[EMu].insert(trigList[EMu].end(),trigEMu,trigEMu+(sizeof(trigEMu)/sizeof(trigEMu[0])));
    trigList[SinglePhoton].insert(trigList[SinglePhoton].end(),trigSinglePhoton,trigSinglePhoton+(sizeof(trigSinglePhoton)/sizeof(trigSinglePhoton[0])));
    switch(trig){
      case DoubleMu:
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1LL<<trigList.at(DoubleMu).at(i))) return TrigHltDiMu_prescale->at(trigList.at(DoubleMu).at(i));
        break;
      case SingleMu:
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1LL<<trigList.at(DoubleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1LL<<trigList.at(SingleMu).at(i))) return TrigHltMu_prescale->at(trigList.at(SingleMu).at(i));
        break;
      case DoubleE: //Includes also HighPtE
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1LL<<trigList.at(DoubleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1LL<<trigList.at(SingleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1LL<<trigList.at(DoubleE).at(i))) return TrigHltDiEl_prescale->at(trigList.at(DoubleE).at(i));
        for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1LL<<trigList.at(HighPtE).at(i))) return TrigHltDiEl_prescale->at(trigList.at(HighPtE).at(i)); //Accepted also for HighPtE
        break;
      case SingleE:
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1LL<<trigList.at(DoubleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1LL<<trigList.at(SingleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1LL<<trigList.at(DoubleE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1LL<<trigList.at(HighPtE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1LL<<trigList.at(SingleE).at(i))) return TrigHltEl_prescale->at(trigList.at(SingleE).at(i));
        break;
      case EMu:
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1LL<<trigList.at(DoubleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1LL<<trigList.at(SingleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1LL<<trigList.at(DoubleE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1LL<<trigList.at(HighPtE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1LL<<trigList.at(SingleE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[EMu].size() ; i++)  if(TrigHltElMu & (1LL<<trigList.at(EMu).at(i))) return TrigHltElMu_prescale->at(trigList.at(EMu).at(i));
        break;
      case SinglePhoton:
        for(unsigned int i = 0 ; i < trigList[SinglePhoton].size() ; i++){
          if(TrigHltPhot & (1LL<<trigList.at(SinglePhoton).at(i))){
            if(selectedPhotonPt){
              if((selectedPhotonPt >= lowThresholdPt[i]) && (selectedPhotonPt <= highThresholdPt[i]+10)) return TrigHltPhot_prescale->at(trigList.at(SinglePhoton).at(i));
            }
            else return TrigHltPhot_prescale->at(trigList.at(SinglePhoton).at(i));
          }
        }
        break;
      case MC_DiLepton://In this case (used for MC ), take if any trigger passed
        for(unsigned int i = 0 ; i < trigList[EMu].size() ; i++)  if(TrigHltElMu & (1LL<<trigList.at(EMu).at(i))) return TrigHltElMu_prescale->at(trigList.at(EMu).at(i));
        for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1LL<<trigList.at(SingleE).at(i))) return TrigHltEl_prescale->at(trigList.at(SingleE).at(i));
        for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1LL<<trigList.at(HighPtE).at(i))) return TrigHltDiEl_prescale->at(trigList.at(HighPtE).at(i));
        for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1LL<<trigList.at(DoubleE).at(i))) return TrigHltDiEl_prescale->at(trigList.at(DoubleE).at(i));
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1LL<<trigList.at(SingleMu).at(i))) return TrigHltMu_prescale->at(trigList.at(SingleMu).at(i));
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1LL<<trigList.at(DoubleMu).at(i))) return TrigHltDiMu_prescale->at(trigList.at(DoubleMu).at(i));
        break;
      case MC_Photon://In this case (used for MC ), take if any trigger passed
        for(unsigned int i = 0 ; i < trigList[SinglePhoton].size() ; i++){
          if(TrigHltPhot & (1LL<<trigList.at(SinglePhoton).at(i))){
            if(selectedPhotonPt){
              if((selectedPhotonPt >= lowThresholdPt[i]) && (selectedPhotonPt <= highThresholdPt[i]+10)) return TrigHltPhot_prescale->at(trigList.at(SinglePhoton).at(i));
            }
            else return TrigHltPhot_prescale->at(trigList.at(SinglePhoton).at(i));
          }
        }
        break;
      default:
        return 0;
    }
    return 0; //If nothing found.
  }

}
