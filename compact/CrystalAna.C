
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBrowser.h"
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include "TH2.h"
#include "TRandom.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Factories.h"
#include "DDG4/Geant4Particle.h"
#include "DDG4/Geant4Data.h"
#include "../src/DualCrystalCalorimeterHit.h"

#include <vector>
#include <algorithm>

const int nchan = 4;
const int ichan[nchan] = {75,73,77,64};
//const int ichan[nchan] = {358,294,422,6};
//const int ichan[nchan] = {326,294,358,6};
//const int ichan[nchan] = {198,166,230,6};
float wavelencut=550;

void crystalana(int num_evtsmax, const char* inputfilename) {


  typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
  typedef std::vector<CalVision::DualCrystalCalorimeterHit*> CalHits;

  // read in libraries that define the classes
  Long_t result;
  char text[1024];
  const char* dd4hep = gSystem->Getenv("DD4hepINSTALL");
  snprintf(text,sizeof(text)," -I%s/include -D__DD4HEP_DDEVE_EXCLUSIVE__ -Wno-shadow -g -O0",dd4hep);
  gSystem->AddIncludePath(text);
  TString fname = "libDDG4IO";
  const char* io_lib = gSystem->FindDynamicLibrary(fname,kTRUE);
  result = gSystem->Load("libDDG4IO");
  result = gSystem->Load("libDDEvePlugins");
  result = gSystem->Load("libDDEvePlugins");
  result = gSystem->Load("libDDDualCrystal");
  result = gSystem->Load("libDDG4Plugins");


  // define histograms
  TH1F *hgenPsize = new TH1F("hgenPsize","number of generator particles",600,0.,40000);
  TH1F *hgenPdgID = new TH1F("hgenpdgID","pdgID of generator particles",600,-200,200);
  TH1F *hcEcalE = new TH1F("hcEcalE","sum crystal ecal energy",100,0.,100.);
  TH1F *hcEcalncer = new TH1F("hcEcalncer","total number of cerenkov",100,0.,10000);
  TH1F *hcEcalncer0 = new TH1F("hcEcalncer0","total number of cerenkov chan 1",100,0.,10000);
  TH1F *hcEcalncer1 = new TH1F("hcEcalncer1","total number of cerenkov chan 2",100,0.,10000);
  TH1F *hcEcalncer2 = new TH1F("hcEcalncer2","total number of cerenkov chan 3",100,0.,10000);
  TH1F *hcEcalncer3 = new TH1F("hcEcalncer3","total number of cerenkov chan 4",100,0.,10000);
  //TH1F *hccerwave = new TH1F("hccerwave","wavelength distribution of cerenkov chan 1",100,0.,10000);
  //TH1F *hcscintwave = new TH1F("hcscintwave","wavelength distribution of scintillation chan 1",100,0.,10000);
  // open data and output file for histograms

  //  const char* inputfilename="/data/users/eno/dd4hep/DD4hep/DDDetectors/compact/testSid.root";
  const char* outputfilename="hist.root";

  // get Tree
  //  TFile *f = new TFile(inputfilename);
  //f->Print();
  GenParts* pgenparts = new GenParts();
  CalHits* pcalhits = new CalHits();
  int num_evt;

  TFile* f = TFile::Open(inputfilename);
  TTree* t = (TTree*)f->Get("EVENT;1");
  t->Print();



  
  // loop over events in the gen loop
  TBranch* b_mc = t->GetBranch("MCParticles");
  int ihaha = b_mc->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<" doing "<<b_mc->GetName()<<std::endl;
  std::cout<<"num_evt gen loop is "<<num_evt<<std::endl;
  
  
  if(num_evt>0) {
    GenParts* gens = new GenParts();
    b_mc->SetAddress(&gens);
    for(int ievt=0;ievt<num_evt; ++ievt) {
      std::cout<<"event number is "<<ievt<<std::endl;
      int nbyte = b_mc->GetEntry(ievt);
      if( nbyte>0) {
	std::cout<<" Gen parts "<<nbyte<<" bytes "<<gens->size() <<std::endl;
      }
      hgenPsize->Fill(gens->size());
      for(size_t i=0;i<gens->size(); ++i) {
        dd4hep::sim::Geant4Particle* agen =gens->at(i);
        hgenPdgID->Fill(agen->pdgID);
      }
    }
  }
  


std::cout<<std::endl;
  

  // loop over events in the ecal loop
  TBranch* b_ecal = t->GetBranch("DRCNoSegment");
  ihaha=b_ecal->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<" doing "<<b_ecal->GetName()<<std::endl;
  std::cout<<"num_evt ecal hit loop is "<<num_evt<<std::endl;
  
  if(num_evt>0) {
    CalHits* ecalhits = new CalHits();
    b_ecal->SetAddress(&ecalhits);

      float esum_avg=0.;
      float esumchan_avg[nchan]={0.,0.,0.,0.};
      float ncerchan_avg[nchan]={0.,0.,0.,0.};
      float nscintchan_avg[nchan]={0.,0.,0.,0.};
      float ncercutchan_avg[nchan]={0.,0.,0.,0.}; //Cerenkov Photons above cutoff
      float nscintcutchan_avg[nchan]={0.,0.,0.,0.}; //Scintillation photons below cutoff
      
      float ncertot_avg=0;
      float nscinttot_avg=0;
      std::vector<float> ncerwavechan_avg[nchan];
      std::vector<float> nscintwavechan_avg[nchan];
      std::vector<float> ncerwavecutchan_avg[nchan]; //Cerenkov Photons above cutoff
      std::vector<float> nscintwavecutchan_avg[nchan]; //Scintillation photons below cutoff
      std::vector<float> number_of_bins_avg[nchan];

    for(int ievt=0;ievt<num_evt; ++ievt) {
      std::cout<<"event number is "<<ievt<<std::endl;
      int nbyte = b_ecal->GetEntry(ievt);
      if( nbyte>0) {
      std::cout<<" Ecal Hits "<<nbyte<<" bytes "<<ecalhits->size() <<std::endl;
      }
      float esum=0.;
      float esumchan[nchan]={0.,0.,0.,0.};
      int ncerchan[nchan]={0,0,0,0};
      int nscintchan[nchan]={0,0,0,0};
      
      //int ncercutchan[nchan]={0,0,0,0}; //Cerenkov Photons above cutoff
      //int nscintcutchan[nchan]={0,0,0,0}; //Scintillation photons below cutoff

      int ncertot=0;
      int nscinttot=0;
      std::vector<int> ncerwavechan[nchan];
      std::vector<int> nscintwavechan[nchan];
      std::vector<int> ncerwavecutchan[nchan]; //Cerenkov Photons above cutoff
      std::vector<int> nscintwavecutchan[nchan]; //Scintillation photons below cutoff      
      std::vector<int> number_of_bins[nchan];
      for(size_t i=0;i<ecalhits->size(); ++i) 
      {
	CalVision::DualCrystalCalorimeterHit* aecalhit =ecalhits->at(i);
	//	std::cout<<"       "<<i<<" energy "<<aecalhit->energyDeposit<<std::endl;
	esum+=aecalhit->energyDeposit;
	ncertot+=aecalhit->ncerenkov;
	nscinttot+=aecalhit->nscintillator;
	//std::cout<<" hit channel is "<< aecalhit->cellID<<" "<<aecalhit->energyDeposit<<" "<<aecalhit->ncerenkov<<" "<<aecalhit->nscintillator<<std::endl;
	int ijchan=aecalhit->nbin;
 	float binsize=(aecalhit->wavelenmax-aecalhit->wavelenmin)/ijchan; //Copied from SDAction
  int bincut=(wavelencut-aecalhit->wavelenmin)/binsize;
  //std::cout<<bincut<<std::endl;
	/*for (int j=0;j<ijchan;j++) {
	  std::cout<<"  ncerwave["<<j<<"]="<<(aecalhit->ncerwave)[j]<<std::endl;
	  std::cout<<"  nscintwave["<<j<<"]="<<(aecalhit->nscintwave)[j]<<std::endl;
	}*/

      // there is a better way to do this
	int jchan=aecalhit->cellID;
	int kchan=-1;
	for( int i=0;i<nchan;i++ ) {
	  if(ichan[i]==jchan) kchan=i; //comparing with the entries of ichan which is populated with the cellIDs 
	}
	if(kchan==-1)
  {
	  //std::cout<<"unknown hit channel is "<< aecalhit->cellID<<std::endl;
	} 
  else
  {
	  esumchan[kchan]+=aecalhit->energyDeposit;
	  ncerchan[kchan]+=aecalhit->ncerenkov;
	  nscintchan[kchan]+=aecalhit->nscintillator;
    for (int j=0; j<ijchan;j++)
    {
     ncerwavechan[kchan].push_back((aecalhit->ncerwave)[j]);
     nscintwavechan[kchan].push_back((aecalhit->nscintwave)[j]);
     number_of_bins[kchan].push_back(j*binsize + aecalhit->wavelenmin);
     if(j<bincut && j>=0)
     {
      ncerwavecutchan[kchan].push_back(0); //'Empty' below the cutoff
      nscintwavecutchan[kchan].push_back((aecalhit->nscintwave)[j]);
     }
     else if (j>=bincut && j<ijchan)
     {
      ncerwavecutchan[kchan].push_back((aecalhit->ncerwave)[j]);
      nscintwavecutchan[kchan].push_back(0); //'Empty' above the cutoff 
     }
    }    
    std::cout<<"Pushed back"<< " "<< ijchan<<std::endl;
	}


      }  // end loop over hits

      esum_avg+=esum;
      ncertot_avg+=ncertot;
      nscinttot_avg+=nscinttot;
      
      for(int k=0;k<nchan;k++)
      {
       esumchan_avg[k]+=esumchan[k];
       ncerchan_avg[k]+=ncerchan[k];
       nscintchan_avg[k]+=nscintchan[k];
       std::cout << ncerwavechan[k].size() << std::endl;
       for(int j=0;j<ncerwavechan[k].size();j++) //Right now all sizes are the same = binsize
       {
        if(ievt==0)
        {
         ncerwavechan_avg[k].push_back(ncerwavechan[k].at(j));
         nscintwavechan_avg[k].push_back(nscintwavechan[k].at(j));
         std::cout<<"Pushed back successfully"<<std::endl;

         ncerwavecutchan_avg[k].push_back(ncerwavecutchan[k].at(j));
         nscintwavecutchan_avg[k].push_back(nscintwavecutchan[k].at(j));
         number_of_bins_avg[k].push_back(number_of_bins[k].at(j));
        }
        else if (ievt>0 &&ievt<num_evt) //i.e. the avg vector already has this size because this condition is necessarily occurring after the 'if' one above
        {
         ncerwavechan_avg[k].at(j) += ncerwavechan[k].at(j);
         nscintwavechan_avg[k].at(j) += nscintwavechan[k].at(j);
         ncerwavecutchan_avg[k].at(j) += (ncerwavecutchan[k].at(j));
         nscintwavecutchan_avg[k].at(j) += (nscintwavecutchan[k].at(j));
        
         //Don't change the number_of_bins_avg here
        }
       }

      }
       
      hcEcalE->Fill(esum/1000.);
      hcEcalncer->Fill(ncertot);
      hcEcalncer0->Fill(ncerchan[0]);
      hcEcalncer1->Fill(ncerchan[1]);
      hcEcalncer2->Fill(ncerchan[2]);
      hcEcalncer3->Fill(ncerchan[3]);

      /*std::cout<<" total energy deposit "<<esum<<std::endl;
      float check=0.;
      for( int i=0;i<nchan;i++) {
	std::cout<<"esum ["<<ichan[i]<<"]="<<esumchan[i]<<std::endl;
	check+=esumchan[i];
      }
      std::cout<<" check total energy deposit "<<check<<std::endl;

      std::cout<<" total number of cherenkov is "<<ncertot<<std::endl;
      check=0;
      for( int i=0;i<nchan;i++) {
	std::cout<<"ncerenkov ["<<ichan[i]<<"]="<<ncerchan[i]<<std::endl;
	check+=ncerchan[i];
      }
      std::cout<<" check ncerenkov "<<check<<std::endl;


      std::cout<<" total number of scintillator is "<<nscinttot<<std::endl;
      check=0;
      for( int i=0;i<nchan;i++) {
	std::cout<<"nscintillator ["<<ichan[i]<<"]="<<nscintchan[i]<<std::endl;
	check+=nscintchan[i];
      }
      std::cout<<" check nscintillator "<<check<<std::endl;
      std::cout<<std::endl;*/

    }  //end loop over events
    
    for( int i=0;i<nchan;i++)
    {
     esumchan_avg[i] = esumchan_avg[i]/num_evt;
     nscintchan_avg[i] = nscintchan_avg[i]/num_evt;
     ncerchan_avg[i] = ncerchan_avg[i]/num_evt;
     for(int j=0;j<ncerwavechan_avg[i].size();j++)
     {
      ncerwavechan_avg[i].at(j) = ncerwavechan_avg[i].at(j)/num_evt;
      nscintwavechan_avg[i].at(j) = nscintwavechan_avg[i].at(j)/num_evt;
      ncerwavecutchan_avg[i].at(j) = ncerwavecutchan_avg[i].at(j)/num_evt;
      nscintwavecutchan_avg[i].at(j) = nscintwavecutchan_avg[i].at(j)/num_evt;
      ncercutchan_avg[i]+= ncerwavecutchan_avg[i].at(j);
      nscintcutchan_avg[i]+= nscintwavecutchan_avg[i].at(j);
     }
     std::cout<<"esum_avg ["<<ichan[i]<<"]="<<esumchan_avg[i]<<std::endl;
     std::cout<<"nscintillator_avg ["<<ichan[i]<<"]="<<nscintchan_avg[i]<<std::endl;
     std::cout<<"ncerenkov_avg ["<<ichan[i]<<"]="<<ncerchan_avg[i]<<std::endl;
     std::cout<<"nscintillator_avg below cutoff ["<<ichan[i]<<"]="<<nscintcutchan_avg[i]<<std::endl;
     std::cout<<"ncerenkov_avg above cutoff ["<<ichan[i]<<"]="<<ncercutchan_avg[i]<<std::endl;     
    }
    
    float norm_cer = *max_element(ncerwavechan_avg[0].begin(), ncerwavechan_avg[0].end());
    float norm_scint = *max_element(nscintwavechan_avg[0].begin(), nscintwavechan_avg[0].end());   
    float norm_cer_det = *max_element(ncerwavechan_avg[2].begin(), ncerwavechan_avg[2].end());
    float norm_scint_det = *max_element(nscintwavechan_avg[2].begin(), nscintwavechan_avg[2].end());   
    
    std::cout<< "max cer is " << " " << norm_cer << endl;
    std::cout<< "max scint is " << " " << norm_scint << endl;
    std::cout<< "max cer det is " << " " << norm_cer_det << endl;
    std::cout<< "max scint det is " << " " << norm_scint_det << endl;
        
    /*for( int i=0;i<nchan;i++)
    {*/
     for(int j=0;j<ncerwavechan_avg[0].size();j++)
     {
      ncerwavechan_avg[0].at(j) = ncerwavechan_avg[0].at(j)/norm_cer;
      nscintwavechan_avg[0].at(j) = nscintwavechan_avg[0].at(j)/norm_scint;
      ncerwavecutchan_avg[0].at(j) = ncerwavecutchan_avg[0].at(j)/norm_cer;
      nscintwavecutchan_avg[0].at(j) = nscintwavecutchan_avg[0].at(j)/norm_scint;
     }
    //}

     for(int j=0;j<ncerwavechan_avg[2].size();j++)
     {
      ncerwavechan_avg[2].at(j) = ncerwavechan_avg[2].at(j)/norm_cer_det/2;
      nscintwavechan_avg[2].at(j) = nscintwavechan_avg[2].at(j)/norm_scint_det/2;
      ncerwavecutchan_avg[2].at(j) = ncerwavecutchan_avg[2].at(j)/norm_cer_det/2;
      nscintwavecutchan_avg[2].at(j) = nscintwavecutchan_avg[2].at(j)/norm_scint_det/2;
     }
     
    TGraph* gcerwave_gen = new TGraph(ncerwavechan_avg[0].size(),&number_of_bins_avg[0][0],&ncerwavechan_avg[0][0]);
    gcerwave_gen->SetTitle("Cerwave_gen");
    TGraph* gcerwave_det = new TGraph(ncerwavechan_avg[1].size(),&number_of_bins_avg[2][0],&ncerwavechan_avg[2][0]);
    gcerwave_det->SetTitle("Cerwave_det");   
     
    TGraph* gscintwave_gen = new TGraph(nscintwavechan_avg[0].size(),&number_of_bins_avg[0][0],&nscintwavechan_avg[0][0]);
    gscintwave_gen->SetTitle("Scintwave_gen");
    TGraph* gscintwave_det = new TGraph(nscintwavechan_avg[1].size(),&number_of_bins_avg[2][0],&nscintwavechan_avg[2][0]);
    gscintwave_det->SetTitle("Scintwave_det");
    
    /*TGraph* gcerwave_gen = new TGraph(ncerwavecutchan_avg[0].size(),&number_of_bins_avg[0][0],&ncerwavecutchan_avg[0][0]);
    gcerwave_gen->SetTitle("Cerwave_gen");
    TGraph* gcerwave_det = new TGraph(ncerwavecutchan_avg[1].size(),&number_of_bins_avg[1][0],&ncerwavecutchan_avg[1][0]);
    gcerwave_det->SetTitle("Cerwave_det");*/   
     
    /*TGraph* gscintwave_gen = new TGraph(nscintwavecutchan_avg[0].size(),&number_of_bins_avg[0][0],&nscintwavecutchan_avg[0][0]);
    gscintwave_gen->SetTitle("Scintwave_gen");
    TGraph* gscintwave_det = new TGraph(nscintwavecutchan_avg[1].size(),&number_of_bins_avg[1][0],&nscintwavecutchan_avg[1][0]);
    gscintwave_det->SetTitle("Scintwave_det");*/ 
    
    TMultiGraph* cerscint = new TMultiGraph();
    cerscint->SetName("cer_scint");
    TCanvas *c1 = new TCanvas("c1","Graph Draw Options", 200,10,600,400); //represents coordinates of start and end points of canvas
    cerscint -> Add(gcerwave_gen);
    cerscint -> Add(gcerwave_det);    
    cerscint -> Add(gscintwave_gen);
    cerscint -> Add(gscintwave_det);
    gcerwave_gen->SetLineColor(kRed);
    gcerwave_det->SetLineColor(kBlack);
    gscintwave_gen->SetLineColor(kBlue);
    gscintwave_det->SetLineColor(kYellow);

    cerscint -> Draw("AC"); //Draw with axes, curve
    TH1F* hist = cerscint->GetHistogram(); //Done purely for tweaking title, axes, etc.
    hist->GetXaxis()->SetTitle("Wavelength (in nm)");
    hist->GetYaxis()->SetTitle("Normalized counts");
    hist->SetTitle("C and S produced in crystal and detected for right killMedia)");
    c1 -> BuildLegend();
    
  }  // end if no events
    
  


  
 
 
  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hgenPsize->Write();
  hgenPdgID->Write();
  hcEcalE->Write();
  hcEcalncer->Write();
  hcEcalncer0->Write();
  hcEcalncer1->Write();
  hcEcalncer2->Write();
  hcEcalncer3->Write();
  //gcerwave->Write();
  //gscintwave->Write();
  out->Close();

}


