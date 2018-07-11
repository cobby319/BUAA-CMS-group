#include <iostream>
#include <thread>
#include <TString.h>
#include <TChain.h>
#include "Includes/ArgParser.h"
#include "Includes/ntuple.h"


using namespace std;

int main(int argc, char **argv)
{
TString InputFile = "OUTPUTS/DoubleMuon_Run2016/CONFIGS/DoubleMuon_jpsipipi.root";
  TString OutputFile = "theOutputFile.root";
  //--- Parse the arguments -----------------------------------------------------
  if (argc > 1) 
  {
    for (int i = 1; i < argc; ++i) 
    {
      TString currentArg = argv[i];
        //--- possible options ---
      if (currentArg.BeginsWith("--input=")) 
      {
        getArg(currentArg, InputFile);
      }
      else if (currentArg.BeginsWith("--output=")) 
      {
        getArg(currentArg, OutputFile);
      }
        
    }
  }


ntuple mytree(InputFile, OutputFile);
mytree.Loop();
return 0;

}
