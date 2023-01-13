
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
//const int ichan[nchan] = {75,73,77,64};
const int ichan[nchan] = {358,294,422,6};
//const int ichan[nchan] = {326,294,358,6};
//const int ichan[nchan] = {198,166,230,6};
float wavelencut=550;

float UV_sipm_QE_x[23] = {361.161, 364.766, 379.794, 387.614,
                           396.624, 406.226, 411.617, 426.594, 436.769, 455.931, 477.492,
                           496.061, 517.627, 547.583, 573.349, 598.521, 615.299, 649.46,
                           671.039, 705.202, 755.548, 773.531, 798.108};
float UV_sipm_QE_y[23] = {0.770120854, 0.787348933, 0.879304547,
                           0.942520324, 0.982752141, 1, 0.982752141, 0.942520324, 0.890796527,
                           0.816088771, 0.741381015, 0.683901339, 0.620685563, 0.545977807,
                           0.488498131, 0.448266313, 0.413790375, 0.356330478, 0.32759064,
                           0.275866843, 0.201139308, 0.178155349, 0.149415511};
                           
                           
float RGB_sipm_QE_x[29] = {305.28, 318.47, 334.67, 352.06,
                            370.06, 396.44, 416.23, 443.81, 466.6, 477.39, 491.78, 515.17,
                            529.56, 556.53, 582.91, 610.49, 636.26, 663.24, 684.22, 712.39,
                            738.76, 755.55, 774.73, 795.11, 825.68, 850.26, 874.23, 894.61, 900.61};
float RGB_sipm_QE_y[29] = {0.034678173, 0.144499016, 0.271678829,
                            0.427750492, 0.525998688, 0.635839415, 0.705195761, 0.786124754,
                            0.87860651, 0.907518244, 0.936410093, 0.994213676, 1, 0.97687459,
                            0.942196417, 0.90173192, 0.849714661, 0.78033843, 0.734107494, 0.664731264,
                            0.583802271, 0.520232248, 0.485554075, 0.427750492, 0.364160585, 0.289017916,
                            0.225428009, 0.167624426, 0.144499016};
//(Directly copied from Yihui's code ana.C)
//Once we get a proper filter with non trivial wavelength dependence it will probably have to be put here

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

      std::vector<float> filter[nchan];
      std::vector<float> QE[nchan];
      //Defined for all channels including the crystal and the air, of course it is for convenience but does not make much physical sense, so all the elements of those two channels will be set to 1 for multiplying i.e. no action

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
     number_of_bins[kchan].push_back(j*binsize + aecalhit->wavelenmin); //lower wavelength of the j th bin
     /*if(j<bincut && j>=0)
     {
      ncerwavecutchan[kchan].push_back(0); //'Empty' below the cutoff
      nscintwavecutchan[kchan].push_back((aecalhit->nscintwave)[j]);
     }
     else if (j>=bincut && j<ijchan)
     {
      ncerwavecutchan[kchan].push_back((aecalhit->ncerwave)[j]);
      nscintwavecutchan[kchan].push_back(0); //'Empty' above the cutoff 
     }*/
    }    
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

         /*ncerwavecutchan_avg[k].push_back(ncerwavecutchan[k].at(j));
         nscintwavecutchan_avg[k].push_back(nscintwavecutchan[k].at(j));*/
         number_of_bins_avg[k].push_back(number_of_bins[k].at(j));
        }
        else if (ievt>0 &&ievt<num_evt) //i.e. the avg vector already has this size because this condition is necessarily occurring after the 'if' one above
        {
         ncerwavechan_avg[k].at(j) += ncerwavechan[k].at(j);
         nscintwavechan_avg[k].at(j) += nscintwavechan[k].at(j);
         /*ncerwavecutchan_avg[k].at(j) += (ncerwavecutchan[k].at(j));
         nscintwavecutchan_avg[k].at(j) += (nscintwavecutchan[k].at(j));*/
        
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
    
    
    
    for(int j=0;j<ncerwavechan_avg[0].size();j++) //Right now all sizes are the same = binsize, taken the channel 0
    {
     if (number_of_bins_avg[1].at(j) <= wavelencut) //Below 550 nm (UV SiPM side, so this is exclusively for scintillation)
     {
      filter[1].push_back(1.0);
     }      
     else if (number_of_bins_avg[1].at(j) > wavelencut)
     {
      filter[1].push_back(0.0);
     }
     
     if (number_of_bins_avg[2].at(j) >= wavelencut) //Above 550 nm (RGB SiPM side, so this is exclusively for Cerenkov)
     {
      filter[2].push_back(1.0);
     }      
     else if (number_of_bins_avg[1].at(j) < wavelencut)
     {
      filter[2].push_back(0.0);
     }
     
     filter[0].push_back(1.0);
     filter[3].push_back(1.0);
    }
    
      for(int j=0;j<ncerwavechan_avg[0].size();j++) //Right now all sizes are the same = binsize, taken the channel 0
      {
      int ent_UV = 22, ent_RGB = 28;
      if(number_of_bins_avg[1].at(j) >= UV_sipm_QE_x[0] && number_of_bins_avg[1].at(j) < UV_sipm_QE_x[ent_UV]) //i.e. if it lies within the range of sensitivity of the UV SiPM
      //1 less than actual size i.e. 23 because of usual indexing quirks (starts from 0)
      {
       for(int i=0;i<ent_UV;i++)
       {
        if (number_of_bins_avg[1].at(j) >= UV_sipm_QE_x[i] && number_of_bins_avg[1].at(j) < UV_sipm_QE_x[i + 1]) //Checking what 'bin' the wavelength lies in 
        QE[1].push_back(0.43 * ((UV_sipm_QE_y[i + 1] - UV_sipm_QE_y[i]) * (number_of_bins_avg[1].at(j) - UV_sipm_QE_x[i]) / (UV_sipm_QE_x[i + 1] - UV_sipm_QE_x[i]) + UV_sipm_QE_y[i]));
       }
      } 
       else 
       {
        QE[1].push_back(0.0);
       }
       
       
       
      if(number_of_bins_avg[2].at(j) >= RGB_sipm_QE_x[0] && number_of_bins_avg[2].at(j) < RGB_sipm_QE_x[ent_RGB]) //i.e. if it lies within the range of sensitivity of the RGB SiPM
      //1 less than actual size i.e. 29 because of usual indexing quirks (starts from 0)
      {
       for(int i=0;i<ent_RGB;i++)
       {
        if (number_of_bins_avg[2].at(j) >= RGB_sipm_QE_x[i] && number_of_bins_avg[2].at(j) < RGB_sipm_QE_x[i + 1]) //Checking what 'bin' the wavelength lies in
        QE[2].push_back(0.325 * ((RGB_sipm_QE_y[i + 1] - RGB_sipm_QE_y[i]) * (number_of_bins_avg[2].at(j) - RGB_sipm_QE_x[i]) / (RGB_sipm_QE_x[i + 1] - RGB_sipm_QE_x[i]) + RGB_sipm_QE_y[i]));
       }
      } 
       else 
       {
        QE[2].push_back(0.0);
       }
       
       QE[0].push_back(1.0);
       QE[3].push_back(1.0);
      }


    for( int i=0;i<nchan;i++)
    {
     esumchan_avg[i] = esumchan_avg[i]/num_evt;
     nscintchan_avg[i] = nscintchan_avg[i]/num_evt;
     ncerchan_avg[i] = ncerchan_avg[i]/num_evt;
     for(int j=0;j<ncerwavechan_avg[i].size();j++)
     {
      ncerwavechan_avg[i].at(j) = ncerwavechan_avg[i].at(j)/num_evt;
      nscintwavechan_avg[i].at(j) = nscintwavechan_avg[i].at(j)/num_evt;
      /*ncerwavecutchan_avg[i].at(j) = ncerwavecutchan_avg[i].at(j)/num_evt;
      nscintwavecutchan_avg[i].at(j) = nscintwavecutchan_avg[i].at(j)/num_evt;*/
      
      ncerwavecutchan_avg[i].push_back(ncerwavechan_avg[i].at(j)*filter[i].at(j));
      nscintwavecutchan_avg[i].push_back(nscintwavechan_avg[i].at(j)*filter[i].at(j));
      
      ncerwavecutchan_avg[i].at(j) = ncerwavecutchan_avg[i].at(j)*QE[i].at(j);
      nscintwavecutchan_avg[i].at(j) = nscintwavecutchan_avg[i].at(j)*QE[i].at(j);
      
      ncercutchan_avg[i]+= ncerwavecutchan_avg[i].at(j);
      nscintcutchan_avg[i]+= nscintwavecutchan_avg[i].at(j);
     }
     std::cout<<"esum_avg ["<<ichan[i]<<"]="<<esumchan_avg[i]<<std::endl;
     std::cout<<"nscintillator_avg ["<<ichan[i]<<"]="<<nscintchan_avg[i]<<std::endl;
     std::cout<<"ncerenkov_avg ["<<ichan[i]<<"]="<<ncerchan_avg[i]<<std::endl;
     std::cout<<"nscintillator_avg below cutoff with QE["<<ichan[i]<<"]="<<nscintcutchan_avg[i]<<std::endl;
     std::cout<<"ncerenkov_avg above cutoff with QE["<<ichan[i]<<"]="<<ncercutchan_avg[i]<<std::endl;     
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
      ncerwavechan_avg[2].at(j) = ncerwavechan_avg[2].at(j)/norm_cer_det/*/2*/;
      nscintwavechan_avg[1].at(j) = nscintwavechan_avg[1].at(j)/norm_scint_det/*/2*/;
      ncerwavecutchan_avg[2].at(j) = ncerwavecutchan_avg[2].at(j)/norm_cer_det/*/2*/;
      nscintwavecutchan_avg[1].at(j) = nscintwavecutchan_avg[1].at(j)/norm_scint_det/*/2*/;
     }
     
    /*TGraph* gcerwave_gen = new TGraph(ncerwavechan_avg[0].size(),&number_of_bins_avg[0][0],&ncerwavechan_avg[0][0]);
    gcerwave_gen->SetTitle("Cerwave_gen");*/
    TGraph* gcerwave_det = new TGraph(ncerwavechan_avg[2].size(),&number_of_bins_avg[2][0],&ncerwavechan_avg[2][0]);
    gcerwave_det->SetTitle("Cerwave_det");   
     
    /*TGraph* gscintwave_gen = new TGraph(nscintwavechan_avg[0].size(),&number_of_bins_avg[0][0],&nscintwavechan_avg[0][0]);
    gscintwave_gen->SetTitle("Scintwave_gen");*/
    TGraph* gscintwave_det = new TGraph(nscintwavechan_avg[1].size(),&number_of_bins_avg[1][0],&nscintwavechan_avg[1][0]);
    gscintwave_det->SetTitle("Scintwave_det");
    
    /*TGraph* gcerwave_gen_cut = new TGraph(ncerwavecutchan_avg[0].size(),&number_of_bins_avg[0][0],&ncerwavecutchan_avg[0][0]);
    gcerwave_gen->SetTitle("Cerwave_gen with cutoff");*/
    TGraph* gcerwave_det_cut = new TGraph(ncerwavecutchan_avg[2].size(),&number_of_bins_avg[2][0],&ncerwavecutchan_avg[2][0]);
    gcerwave_det_cut->SetTitle("Cerwave_det with cutoff and QE");
     
    /*TGraph* gscintwave_gen_cut = new TGraph(nscintwavecutchan_avg[0].size(),&number_of_bins_avg[0][0],&nscintwavecutchan_avg[0][0]);
    gscintwave_gen->SetTitle("Scintwave_gen with cutoff");*/
    TGraph* gscintwave_det_cut = new TGraph(nscintwavecutchan_avg[1].size(),&number_of_bins_avg[1][0],&nscintwavecutchan_avg[1][0]);
    gscintwave_det_cut->SetTitle("Scintwave_det with cutoff and QE"); 
    
    TGraph* QE_Crystal = new TGraph(QE[0].size(),&number_of_bins_avg[0][0],&QE[0][0]);
    QE_Crystal->SetTitle("QE for PbWO4 crystal");
    TGraph* QE_Air_out = new TGraph(QE[3].size(),&number_of_bins_avg[3][0],&QE[3][0]);
    QE_Air_out->SetTitle("QE for Air outside");
    
    TGraph* QE_UV = new TGraph(QE[1].size(),&number_of_bins_avg[1][0],&QE[1][0]);
    QE_UV->SetTitle("QE for Left killMedia (UV)");
    TGraph* QE_RGB = new TGraph(QE[2].size(),&number_of_bins_avg[2][0],&QE[2][0]);
    QE_RGB->SetTitle("QE for Right killMedia (RGB)");  
    
    
    TGraph* filter_Crystal = new TGraph(filter[0].size(),&number_of_bins_avg[0][0],&filter[0][0]);
    filter_Crystal->SetTitle("filter for PbWO4 crystal");
    TGraph* filter_Air_out = new TGraph(filter[3].size(),&number_of_bins_avg[3][0],&filter[3][0]);
    filter_Air_out->SetTitle("filter for Air outside");
    
    TGraph* filter_UV = new TGraph(filter[1].size(),&number_of_bins_avg[1][0],&filter[1][0]);
    filter_UV->SetTitle("filter for Left killMedia (UV)");
    TGraph* filter_RGB = new TGraph(filter[2].size(),&number_of_bins_avg[2][0],&filter[2][0]);
    filter_RGB->SetTitle("filter for Right killMedia (RGB)");  

    TMultiGraph* cerscint = new TMultiGraph();
    cerscint->SetName("cer_scint");
    TCanvas *c1 = new TCanvas("c1","Graph Draw Options", 200,10,600,400); //represents coordinates of start and end points of canvas
    //cerscint -> Add(gcerwave_gen);
    cerscint -> Add(gcerwave_det);    
    //cerscint -> Add(gscintwave_gen);
    cerscint -> Add(gscintwave_det);
    cerscint -> Add(gcerwave_det_cut);
    cerscint -> Add(gscintwave_det_cut);
    //gcerwave_gen->SetLineColor(kRed);
    gcerwave_det->SetLineColor(kBlack);
    //gscintwave_gen->SetLineColor(kBlue);
    gscintwave_det->SetLineColor(kOrange);
    gcerwave_det_cut->SetLineColor(kRed);
    gscintwave_det_cut->SetLineColor(kBlue);
    
    /*cerscint -> Add(QE_Crystal);
    cerscint -> Add(QE_Air_out);
    cerscint -> Add(QE_UV);
    cerscint -> Add(QE_RGB);
    QE_Crystal->SetLineColor(kRed);
    QE_Air_out->SetLineColor(kBlack);
    QE_UV->SetLineColor(kBlue);
    QE_RGB->SetLineColor(kOrange);*/
    
    /*cerscint -> Add(filter_Crystal);
    cerscint -> Add(filter_Air_out);
    cerscint -> Add(filter_UV);
    cerscint -> Add(filter_RGB);
    //filter_Crystal->SetMarkerColor(kRed);
    //filter_Air_out->SetMarkerColor(kBlack);
    //filter_UV->SetLineColor(kBlue);
    filter_UV->SetMarkerColor(kBlue);
    //filter_RGB->SetLineColor(kOrange);
    filter_RGB->SetMarkerColor(kOrange);*/
     
    cerscint -> Draw("AC"); //Draw with axes, curve
    TH1F* hist = cerscint->GetHistogram(); //Done purely for tweaking title, axes, etc.
    hist->GetXaxis()->SetTitle("Wavelength (in nm)");
    hist->GetYaxis()->SetTitle("Normalized counts");
    hist->SetTitle("C (right killMedia) and S (left killMedia) with cutoff and QE)");
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


