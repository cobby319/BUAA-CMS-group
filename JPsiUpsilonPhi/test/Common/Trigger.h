#ifndef trigger_h
#define trigger_h

#include <iostream>
#include <string>
#include <TMath.h>

namespace trigger
{
  enum {DoubleMu, SingleMu, DoubleE, HighPtE, SingleE, EMu, SinglePhoton, MC_DiLepton, MC_Photon}; //List of triggers used for our analysis

  int passTrigger(int trig, ULong64_t TrigHltDiMu, ULong64_t TrigHltMu, ULong64_t TrigHltDiEl, ULong64_t TrigHltEl, ULong64_t TrigHltElMu, ULong64_t TrigHltPhot, std::vector<unsigned int> *TrigHltDiMu_prescale, std::vector<unsigned int> *TrigHltMu_prescale, std::vector<unsigned int> *TrigHltDiEl_prescale, std::vector<unsigned int> *TrigHltEl_prescale, std::vector<unsigned int> *TrigHltElMu_prescale, std::vector<unsigned int> *TrigHltPhot_prescale, double selectedPhotonPt = 0);
}

#endif
