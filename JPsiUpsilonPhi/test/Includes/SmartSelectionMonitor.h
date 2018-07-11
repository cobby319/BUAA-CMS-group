#ifndef smartselectionmonitor_hh
#define smartselectionmonitor_hh

// system include files
#include<iostream>
#include<string>
#include<map>
#include<unordered_map>
#include<algorithm>
#include<vector>
#include<memory>

// user include files
#include "TH1D.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TString.h"
#include "TROOT.h"

namespace std{
  template<> struct hash< TString >{ size_t operator()( const TString& x ) const{ return hash<std::string>()( x.Data() );  }  };
}

class SmartSelectionMonitor {
  
public:

  SmartSelectionMonitor(){}
  ~SmartSelectionMonitor() { }


  //types
  typedef std::unordered_map<TString, std::map<TString, TH1*>* > Monitor_t;


  //short getters
  inline Monitor_t &getAllMonitors() { return allMonitors_; }

  //checks if base Histo Exist
  inline bool hasBaseHisto(TString histo){
    if(allMonitors_.find(histo) == allMonitors_.end())return false;
    return true;
  }

  //checks if tag Exist for a given histo
  inline bool hasTag(std::map<TString, TH1*>* map, TString tag, bool useBinWidth = false){
    if( map->find(tag) != map->end() )return true;
    if( map->find("all") == map->end() )return false;
  
    TH1* base = (*map)["all"];
    TString allName = base->GetName();
    TString name = allName + "_" + tag;
    TH1* h = (TH1*) base->Clone(name);
    h->SetName(name);
    h->SetTitle(name);
    if(useBinWidth){
      TString yaxisTitle = h->GetYaxis()->GetTitle();
      if(!yaxisTitle.Contains("bin")) yaxisTitle += " / GeV";
      h->GetYaxis()->SetTitle(yaxisTitle);
    }
    h->Reset("ICE");
    h->SetDirectory(gROOT);
    (*map)[tag] = h; 
    //printf("new histo created with name = %30s and tag = %15s: Name=%s\n",allName.Data(), tag.Data(), h->GetName());
    return true;
  }
  
  //get histo
  inline TH1 *getHisto(TString histo,TString tag, bool useBinWidth = false){
    if( !hasBaseHisto(histo) )return NULL;
    std::map<TString, TH1*>* map = allMonitors_[histo];
    if( !hasTag(map, tag, useBinWidth) )return NULL;
    return (*map)[tag];
  }

  //write all histo
  inline void Write(){
    for(Monitor_t::iterator it =allMonitors_.begin(); it!= allMonitors_.end(); it++){
      std::map<TString, TH1*>* map = it->second;
      bool neverFilled = true;

      for(std::map<TString, TH1*>::iterator h =map->begin(); h!= map->end(); h++){
        if(!(h->second)){printf("histo = %30s %15s IS NULL",it->first.Data(), h->first.Data());continue;}
        if(h->second->GetEntries()>0)neverFilled = false;

        if(h->first=="all"){h->second->SetName(h->first+"_"+h->second->GetName());}
        //printf("histo = %30s tag = %15s Name = %s\n",it->first.Data(), h->first.Data(),  h->second->GetName());
        if(h->first!="all") h->second->Write();
      }

      if(neverFilled){printf("SmartSelectionMonitor: histo = '%s' is empty for all categories, you may want to cleanup your project to remove this histogram\n",it->first.Data());}
    } 
  }

  //scale all histo by w
  inline void Scale(double w){
     for(Monitor_t::iterator it =allMonitors_.begin(); it!= allMonitors_.end(); it++){
        std::map<TString, TH1*>* map = it->second;
        for(std::map<TString, TH1*>::iterator h =map->begin(); h!= map->end(); h++){
	  if(!(h->second)){continue;}
          h->second->Scale(w);
        }
     } 
  }


  //takes care of filling an histogram
  bool fillHisto  (TString name, TString tag, double valx, double weight, bool useBinWidth=false);
  bool fillHisto  (TString name, TString tag, double valx, double valy, double weight,  bool useBinWidth=false);
  bool fillProfile(TString name, TString tag, double valx, double valy, double weight);

  bool fillHisto(TString name, std::vector<TString> tags, double valx, double weight,  bool useBinWidth=false);
  bool fillHisto(TString name, std::vector<TString> tags, double valx, double valy, double weight,  bool useBinWidth=false);
  bool fillProfile(TString name, std::vector<TString> tags, double valx, double valy, double weight);

  bool fillHisto(TString name, std::vector<TString> tags, double valx, std::vector<double> weights,  bool useBinWidth=false);
  bool fillHisto(TString name, std::vector<TString> tags, double valx, double valy, std::vector<double> weights,  bool useBinWidth=false);
  bool fillProfile(TString name, std::vector<TString> tags, double valx, double valy, std::vector<double> weights);


   //short inits the monitor plots for a new step
  void initMonitorForStep(TString tag);
  
  //short add new histogram
  TH1 * addHistogram(TH1 *h, TString tag);
  TH1 * addHistogram(TH1 *h);
  
public: //I know, it's bad. But I didn't find any other way.

  //all the selection step monitors
  Monitor_t allMonitors_;
};

#endif
