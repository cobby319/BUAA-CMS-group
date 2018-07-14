#include <iostream>
#include <thread>
#include <TString.h>
#include <TChain.h>
#include "Includes/ArgParser.h"
#include "Includes/ntuple.h"


using namespace std;

int main(int argc, char **argv)
{
TString InputFile = "/pnfs/iihe/cms/store/user/hanwen/Charmonium/crab_Charmonium_Run2016H-03Feb2017_ver3-v1/180713_140833/0000/Charmonium_Run2016_84.root";
  TString OutputFile = "theOutputFile.root";
  TString catalogInputFile = "";
  int skipFile =0;
  int maxFile = 0;
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
      else if (currentArg.BeginsWith("catalogInputFile=")) {
        getArg(currentArg, catalogInputFile);
      }
      else if (currentArg.BeginsWith("skip-files=")) {
        getArg(currentArg, skipFile);
      }
      else if (currentArg.BeginsWith("max-files=")) {
        getArg(currentArg, maxFile);
      }
      else if (currentArg.BeginsWith("--output=")) 
      {
        getArg(currentArg, OutputFile);
      }
        
    }
  }


ntuple mytree(catalogInputFile, OutputFile, skipFile, maxFile);
mytree.Loop();
return 0;

}
