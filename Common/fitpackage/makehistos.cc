#include "Sources/makeHistos.C"
#include "Includes/ArgParser.h"
#include <iostream>
#include "TString.h"
int main (int argc, char ** argv) {
   int ptbin = 0;
   int etabin = 0;
   if (argc>1) {
      for (int i=0;i<argc;i++){
         TString currentArg = argv[i];
         if (currentArg.BeginsWith("pt")) getArg(currentArg,ptbin);
         
         else if (currentArg.BeginsWith("eta")) getArg(currentArg,etabin);
         }
      }
   std::cout<<"pt&eta bin =="<<ptbin<<"\t"<<etabin<<std::endl;
   makeHistos("merged_DY.root",ptbin,etabin);
   return 0;
   } 
