
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
#include <TPaveText.h>
#include "TRandom.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Factories.h"
#include "DDG4/Geant4Particle.h"
#include "DDG4/Geant4Data.h"
#include "../src/DualCrystalCalorimeterHit.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>

// information about hit channel IT numbers
const int nchan = 4;
//const int ichan[nchan] = {75,73,77,64};
const int ichan[nchan] = {147,145,149,128};
//const int ichan[nchan] = {149,145,153,128};
//const int ichan[nchan] = {358,294,422,6};
//const int ichan[nchan] = {326,294,358,6};
//const int ichan[nchan] = {198,166,230,6};

/*double UV_sipm_QE_x[23] = {361.161, 364.766, 379.794, 387.614, 396.624, 406.226, 411.617, 426.594, 436.769, 455.931, 477.492, 496.061, 517.627, 547.583, 573.349, 598.521, 615.299, 649.46, 671.039, 705.202, 755.548, 773.531, 798.108};

double UV_sipm_QE_y[23] = {0.770120854, 0.787348933, 0.879304547, 0.942520324, 0.982752141, 1, 0.982752141, 0.942520324, 0.890796527, 0.816088771, 0.741381015, 0.683901339, 0.620685563, 0.545977807, 0.488498131, 0.448266313, 0.413790375, 0.356330478, 0.32759064, 0.275866843, 0.201139308, 0.178155349, 0.149415511};
   
double RGB_sipm_QE_x[29] = {305.28, 318.47, 334.67, 352.06, 370.06, 396.44, 416.23, 443.81, 466.6, 477.39, 491.78, 515.17, 529.56, 556.53, 582.91, 610.49, 636.26, 663.24, 684.22, 712.39, 738.76, 755.55, 774.73, 795.11, 825.68, 850.26, 874.23, 894.61, 900.61};

double RGB_sipm_QE_y[29] = {0.034678173, 0.144499016, 0.271678829, 0.427750492, 0.525998688, 0.635839415, 0.705195761, 0.786124754, 0.87860651, 0.907518244, 0.936410093, 0.994213676, 1, 0.97687459, 0.942196417, 0.90173192, 0.849714661, 0.78033843, 0.734107494, 0.664731264, 0.583802271, 0.520232248, 0.485554075, 0.427750492, 0.364160585, 0.289017916, 0.225428009, 0.167624426, 0.144499016};*/
                       
//(Directly copied from Yihui's code ana.C)

double UV_sipm_QE_x[89] = {282.6086957, 288.0434783, 290.2173913, 291.3043478, 292.3913043, 293.4782609, 294.5652174, 295.6521739, 296.7391304, 297.826087, 298.9130435, 301.0869565, 303.2608696, 304.3478261, 307.6086957, 309.7826087, 314.1304348, 317.3913043, 321.7391304, 329.3478261, 335.8695652, 343.4782609, 350., 358.6956522, 367.3913043, 373.9130435, 378.2608696, 383.6956522, 388.0434783, 394.5652174, 402.173913, 408.6956522, 417.3913043, 428.2608696, 440.2173913, 452.173913, 465.2173913, 475., 482.6086957, 489.1304348, 494.5652174, 500., 505.4347826, 510.8695652, 517.3913043, 522.826087, 529.3478261, 534.7826087, 542.3913043, 548.9130435, 555.4347826, 561.9565217, 568.4782609, 577.173913, 583.6956522, 590.2173913, 596.7391304, 603.2608696, 610.8695652, 617.3913043, 623.9130435, 631.5217391, 638.0434783, 644.5652174, 652.173913, 659.7826087, 667.3913043, 676.0869565, 685.8695652, 694.5652174, 703.2608696, 713.0434783, 721.7391304, 731.5217391, 741.3043478, 753.2608696, 764.1304348, 775., 785.8695652, 797.826087, 808.6956522, 819.5652174, 831.5217391, 844.5652174, 855.4347826, 866.3043478, 877.173913, 889.1304348, 898.9130435};

double UV_sipm_QE_y[89] = {0.02983683, 0.041958042, 0.068997669, 0.106293706, 0.131468531, 0.158508159, 0.185547786, 0.211655012, 0.225641026, 0.240559441, 0.255477855, 0.271328671, 0.287179487, 0.299300699, 0.30955711, 0.321678322, 0.331934732, 0.342191142, 0.353379953, 0.363636364, 0.374825175, 0.385081585, 0.393473193, 0.399067599, 0.407459207, 0.418648019, 0.431701632, 0.443822844, 0.455011655, 0.466200466, 0.476456876, 0.485780886, 0.494172494, 0.502564103, 0.508158508, 0.509090909, 0.506293706, 0.498834499, 0.48951049, 0.479254079, 0.468997669, 0.458741259, 0.446620047, 0.437296037, 0.427039627, 0.416783217, 0.405594406, 0.395337995, 0.384149184, 0.373892774, 0.363636364, 0.352447552, 0.342191142, 0.331002331, 0.320745921, 0.31048951, 0.301165501, 0.290909091, 0.280652681, 0.271328671, 0.261072261, 0.251748252, 0.242424242, 0.233100233, 0.222843823, 0.214452214, 0.204195804, 0.193939394, 0.185547786, 0.176223776, 0.167832168, 0.158508159, 0.151048951, 0.142657343, 0.134265734, 0.125874126, 0.118414918, 0.11002331, 0.102564103, 0.095104895, 0.087645688, 0.082051282, 0.074592075, 0.068065268, 0.061538462, 0.055944056, 0.049417249, 0.044755245, 0.038694639};

double RGB_sipm_QE_x[89] = {282.6086957, 288.0434783, 290.2173913, 291.3043478, 292.3913043, 293.4782609, 294.5652174, 295.6521739, 296.7391304, 297.826087, 298.9130435, 301.0869565, 303.2608696, 304.3478261, 307.6086957, 309.7826087, 314.1304348, 317.3913043, 321.7391304, 329.3478261, 335.8695652, 343.4782609, 350., 358.6956522, 367.3913043, 373.9130435, 378.2608696, 383.6956522, 388.0434783, 394.5652174, 402.173913, 408.6956522, 417.3913043, 428.2608696, 440.2173913, 452.173913, 465.2173913, 475., 482.6086957, 489.1304348, 494.5652174, 500., 505.4347826, 510.8695652, 517.3913043, 522.826087, 529.3478261, 534.7826087, 542.3913043, 548.9130435, 555.4347826, 561.9565217, 568.4782609, 577.173913, 583.6956522, 590.2173913, 596.7391304, 603.2608696, 610.8695652, 617.3913043, 623.9130435, 631.5217391, 638.0434783, 644.5652174, 652.173913, 659.7826087, 667.3913043, 676.0869565, 685.8695652, 694.5652174, 703.2608696, 713.0434783, 721.7391304, 731.5217391, 741.3043478, 753.2608696, 764.1304348, 775., 785.8695652, 797.826087, 808.6956522, 819.5652174, 831.5217391, 844.5652174, 855.4347826, 866.3043478, 877.173913, 889.1304348, 898.9130435};

double RGB_sipm_QE_y[89] = {0.02983683, 0.041958042, 0.068997669, 0.106293706, 0.131468531, 0.158508159, 0.185547786, 0.211655012, 0.225641026, 0.240559441, 0.255477855, 0.271328671, 0.287179487, 0.299300699, 0.30955711, 0.321678322, 0.331934732, 0.342191142, 0.353379953, 0.363636364, 0.374825175, 0.385081585, 0.393473193, 0.399067599, 0.407459207, 0.418648019, 0.431701632, 0.443822844, 0.455011655, 0.466200466, 0.476456876, 0.485780886, 0.494172494, 0.502564103, 0.508158508, 0.509090909, 0.506293706, 0.498834499, 0.48951049, 0.479254079, 0.468997669, 0.458741259, 0.446620047, 0.437296037, 0.427039627, 0.416783217, 0.405594406, 0.395337995, 0.384149184, 0.373892774, 0.363636364, 0.352447552, 0.342191142, 0.331002331, 0.320745921, 0.31048951, 0.301165501, 0.290909091, 0.280652681, 0.271328671, 0.261072261, 0.251748252, 0.242424242, 0.233100233, 0.222843823, 0.214452214, 0.204195804, 0.193939394, 0.185547786, 0.176223776, 0.167832168, 0.158508159, 0.151048951, 0.142657343, 0.134265734, 0.125874126, 0.118414918, 0.11002331, 0.102564103, 0.095104895, 0.087645688, 0.082051282, 0.074592075, 0.068065268, 0.061538462, 0.055944056, 0.049417249, 0.044755245, 0.038694639};

//Once we get a proper filter with non trivial wavelength dependence it will probably have to be put here
double UV_filter_x[92] = {300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, 410.0, 420.0, 430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0, 500.0, 510.0, 520.0, 530.0, 540.0, 550.0, 560.0, 570.0, 580.0, 590.0, 600.0, 610.0, 620.0, 630.0, 640.0, 650.0, 660.0, 670.0, 680.0, 690.0, 700.0, 710.0, 720.0, 730.0, 740.0, 750.0, 760.0, 770.0, 780.0, 790.0, 800.0, 810.0, 820.0, 830.0, 840.0, 850.0, 860.0, 870.0, 880.0, 890.0, 900.0, 910.0, 920.0, 930.0, 940.0, 950.0, 960.0, 970.0, 980.0, 990.0, 1000.0, 1010.0, 1020.0, 1030.0, 1040.0, 1050.0, 1060.0, 1070.0, 1080.0, 1090.0, 1100.0, 1110.0, 1120.0, 1130.0, 1140.0, 1150.0, 1160.0, 1170.0, 1180.0, 1190.0, 1200.0, 1250.0};

double UV_filter_y[92] = {0.88, 0.91, 0.9200, 0.9260, 0.9210, 0.9060, 0.8660, 0.7560, 0.4820, 0.1200, 0.003999, 0.00001047, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0003631, 0.008128, 0.05000, 0.1380, 0.2220, 0.2780, 0.2830, 0.2400, 0.1850, 0.1360, 0.09700, 0.06700, 0.04700, 0.03300, 0.02400, 0.01800, 0.01382, 0.01088, 0.009121, 0.007965, 0.007014, 0.006255, 0.005892, 0.005834, 0.005948, 0.006114, 0.006368, 0.007000, 0.008090, 0.009390, 0.01100, 0.01300, 0.01520, 0.01825, 0.02132, 0.02498, 0.02809, 0.03263, 0.03673, 0.04100, 0.04375, 0.04599, 0.04800, 0.04828, 0.04897, 0.04608, 0.04300, 0.03948, 0.03666, 0.03325, 0.02970, 0.02691, 0.02354, 0.02064, 0.01775, 0.008800};               

double GG_filter_x[92] = {300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, 410.0, 420.0, 430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0, 500.0, 510.0, 520.0, 530.0, 540.0, 550.0, 560.0, 570.0, 580.0, 590.0, 600.0, 610.0, 620.0, 630.0, 640.0, 650.0, 660.0, 670.0, 680.0, 690.0, 700.0, 710.0, 720.0, 730.0, 740.0, 750.0, 760.0, 770.0, 780.0, 790.0, 800.0, 810.0, 820.0, 830.0, 840.0, 850.0, 860.0, 870.0, 880.0, 890.0, 900.0, 910.0, 920.0, 930.0, 940.0, 950.0, 960.0, 970.0, 980.0, 990.0, 1000.0, 1010.0, 1020.0, 1030.0, 1040.0, 1050.0, 1060.0, 1070.0, 1080.0, 1090.0, 1100.0, 1110.0, 1120.0, 1130.0, 1140.0, 1150.0, 1160.0, 1170.0, 1180.0, 1190.0, 1200.0, 1250.0};

/*double GG_filter_y[92] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00007048, 0.06112, 0.4792, 0.7703, 0.8551, 0.8849, 0.9032, 0.9176, 0.9300, 0.9409, 0.9500, 0.9578, 0.9645, 0.9703, 0.9747, 0.9781, 0.9809, 0.9829, 0.9843, 0.9852, 0.9857, 0.9858, 0.9856, 0.9852, 0.9847, 0.9841, 0.9834, 0.9825, 0.9814, 0.9803, 0.9792, 0.9780, 0.9767, 0.9754, 0.9740, 0.9724, 0.9706, 0.9691, 0.9676, 0.9662, 0.9649, 0.9635, 0.9622, 0.9609, 0.9596, 0.9585, 0.9574, 0.9563, 0.9552, 0.9541, 0.9531, 0.9523, 0.9517, 0.9512, 0.9507, 0.9503, 0.9499, 0.9495, 0.9491, 0.9488, 0.9484, 0.9481, 0.9478, 0.9475, 0.9472, 0.9470, 0.9468, 0.9467, 0.9466, 0.9465, 0.9465, 0.9465, 0.9465, 0.9465, 0.9466, 0.9467, 0.9468, 0.9469, 0.9470, 0.9472, 0.9474, 0.9484};*/ // This one is for GG420

double GG_filter_y[92] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01043, 0.2606, 0.7129, 0.9022, 0.9523, 0.9691, 0.9773, 0.9819, 0.9850, 0.9866, 0.9877, 0.9882, 0.9881, 0.9876, 0.9872, 0.9864, 0.9857, 0.9850, 0.9842, 0.9833, 0.9822, 0.9811, 0.9799, 0.9783, 0.9769, 0.9754, 0.9739, 0.9722, 0.9704, 0.9686, 0.9667, 0.9649, 0.9632, 0.9616, 0.9600, 0.9583, 0.9568, 0.9554, 0.9541, 0.9529, 0.9516, 0.9504, 0.9492, 0.9482, 0.9474, 0.9468, 0.9461, 0.9454, 0.9449, 0.9443, 0.9438, 0.9434, 0.9429, 0.9424, 0.9420, 0.9416, 0.9412, 0.9409, 0.9407, 0.9406, 0.9405, 0.9404, 0.9404, 0.9404, 0.9404, 0.9404, 0.9404, 0.9405, 0.9406, 0.9407, 0.9408, 0.9409, 0.9411, 0.9413, 0.9415, 0.9430}; //This one is for GG475 which is the one being used

void crystalana(int num_evtsmax, const char* inputfilename) {


  typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
  typedef std::vector<CalVision::DualCrystalCalorimeterHit*> CalHits;
  typedef dd4hep::sim::Geant4HitData::MonteCarloContrib Contribution;
  typedef std::vector<dd4hep::sim::Geant4HitData::MonteCarloContrib> Contributions;


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
  result = gSystem->Load("libSingleDualCrystal");
  result = gSystem->Load("libDDG4Plugins");


  // define histograms
  
  //gen particles
  TH1F *hgenPsize = new TH1F("hgenPsize","number of generator particles",600,0.,40000);
  TH1F *hgenPdgID = new TH1F("hgenpdgID","pdgID of generator particles",600,-200,200);
  
  // calorimeter infor
  TH1F *hchan = new TH1F("hchan","channel ID number",1028,0.,1028);
  TH1F *hcEcalE = new TH1F("hcEcalE","sum crystal ecal energy",100,0.,100.);
  TH1F *hcEcalncer = new TH1F("hcEcalncer","total number of cerenkov",100,0.,10000);
  TH1F *hcEcalncer0 = new TH1F("hcEcalncer0","total number of cerenkov chan 1",100,0.,10000);
  TH1F *hcEcalncer1 = new TH1F("hcEcalncer1","total number of cerenkov chan 2",100,0.,10000);
  TH1F *hcEcalncer2 = new TH1F("hcEcalncer2","total number of cerenkov chan 3",100,0.,10000);
  TH1F *hcEcalncer3 = new TH1F("hcEcalncer3","total number of cerenkov chan 4",100,0.,10000);
  //TH1F *hccerwave = new TH1F("hccerwave","wavelength distribution of cerenkov chan 1",100,0.,10000);
  //TH1F *hcscintwave = new TH1F("hcscintwave","wavelength distribution of scintillation chan 1",100,0.,10000);
  
  TH1F *hcEcalE0 = new TH1F("hcEcalE0","energy chan 0",100,0.,10000);
  TH1F *hcEcalE1 = new TH1F("hcEcalE1","energy chan 1",100,0.,10000);
  TH1F *hcEcalE2 = new TH1F("hcEcalE2","energy chan 2",100,0.,10000);
  TH1F *hcEcalE3 = new TH1F("hcEcalE3","energy chan 3",100,0.,10000);

  // open data and output file for histograms

  //  const char* inputfilename="/data/users/eno/dd4hep/DD4hep/DDDetectors/compact/testSid.root";
  const char* outputfilename="hist.root";

  // get Tree
  //  TFile *f = new TFile(inputfilename);
  //f->Print();
  GenParts* pgenparts = new GenParts();
  CalHits* pcalhits = new CalHits();
  int num_evt,nbyte;

  TFile* f = TFile::Open(inputfilename);
  TTree* t = (TTree*)f->Get("EVENT;1");
  t->Print();



  
  // loop over events
  TBranch* b_mc = t->GetBranch("MCParticles");
  TBranch* b_ecal = t->GetBranch("DRCNoSegment");
  int ihaha = b_mc->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<" doing "<<b_mc->GetName()<<std::endl;
  std::cout<<"num_evt gen loop is "<<num_evt<<std::endl;
  

std::cout<<std::endl;
  
  
  if(num_evt>0) {
  
    // find branches
    GenParts* gens = new GenParts();
    b_mc->SetAddress(&gens);
    CalHits* ecalhits = new CalHits();
    b_ecal->SetAddress(&ecalhits);

      float esum_avg=0.;
      float esumchan_avg[nchan]={0.,0.,0.,0.};
      float ncerchan_avg[nchan]={0.,0.,0.,0.};
      float nscintchan_avg[nchan]={0.,0.,0.,0.};
      
      float ncertot_avg=0;
      float nscinttot_avg=0;
      std::vector<float> ncerwavechan_avg[nchan];
      std::vector<float> nscintwavechan_avg[nchan];
      std::vector<float> ncerwavecutchan_avg[nchan]; //Cerenkov Photons above cutoff
      std::vector<float> nscintwavecutchan_avg[nchan]; //Scintillation photons below cutoff
      std::vector<float> number_of_bins_avg[nchan];
      float a = 2.0;
      int kink_cut;
      /*std::vector<float> ncertimechan_avg[nchan];
      std::vector<float> nscinttimechan_avg[nchan];   
      std::vector<float> number_of_bins_time_avg[nchan];*/

      std::vector<float> filter[nchan];
      std::vector<float> QE[nchan];
      //Defined for all channels including the crystal and the air, of course it is for convenience but does not make much physical sense, so all the elements of those two channels will be set to 1 for multiplying i.e. no action

    int SCEPRINT2=10;
    
    for(int ievt=0;ievt<num_evt; ++ievt) {
      std::cout<<"event number is "<<ievt<<std::endl;
      //gen particles
      nbyte = b_mc->GetEntry(ievt);
      
      if( nbyte>0) {
	if(ievt<SCEPRINT2) std::cout<<gens->size()<<" Gen particles "<<std::endl;
      }
      float mainee=-1;
      hgenPsize->Fill(gens->size());
      for(size_t i=0;i<gens->size(); ++i) {
        dd4hep::sim::Geant4Particle* agen =gens->at(i);
        hgenPdgID->Fill(agen->pdgID);

	float px=agen->psx;
	float py=agen->psy;
	float pz=agen->psz;
	float mass=agen->mass;
	float ee=sqrt(mass*mass+px*px+py*py+pz*pz);

	if(i==0&&ievt<SCEPRINT2) std::cout<<" first particle pid energy are "<<agen->pdgID<<" "<<ee<<std::endl;

	if(ee>mainee) {
	  mainee=ee;
	  if(ievt<SCEPRINT2) std::cout<<"highest e  gen now has pid "<<agen->pdgID<<" energy "<<ee<<std::endl;
	}
	float vsy=agen->vsy;
	/*if(vsy>20)*/ std::cout<<" pid "<<agen->pdgID<<" "
			      <<"vs ("<<agen->vsx<<","<<agen->vsy<<","<<agen->vsz<<") "
			      <<"ve ("<<agen->vex<<","<<agen->vey<<","<<agen->vez<<") "
			      <<"ps ("<<agen->psx<<","<<agen->psy<<","<<agen->psz<<") "
			      <<"pe ("<<agen->pex<<","<<agen->pey<<","<<agen->pez<<") "
			      <<std::endl;



      }



    // ECAL hits  
    // there are hits in the crystal and also the photodetectors "kill media"
    // in the crystal, photons created in the crystal are counted and their wavelengths stored
    // in the photodetector, photons that enter are counted, wavelength stored, and then they are killed


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
      int SCEPRINT=10;
      std::vector<int> ncerwavechan[nchan];
      std::vector<int> nscintwavechan[nchan];
      std::vector<int> ncerwavecutchan[nchan]; //Cerenkov Photons above cutoff
      std::vector<int> nscintwavecutchan[nchan]; //Scintillation photons below cutoff      
      std::vector<float> number_of_bins[nchan];
      
      /*std::vector<int> ncertimechan[nchan];
      std::vector<int> nscinttimechan[nchan];   
      std::vector<float> number_of_bins_time[nchan];*/

      for(size_t i=0;i<ecalhits->size(); ++i) 
      {
	CalVision::DualCrystalCalorimeterHit* aecalhit =ecalhits->at(i);
	//	std::cout<<"       "<<i<<" energy "<<aecalhit->energyDeposit<<std::endl;
	esum+=aecalhit->energyDeposit;
	ncertot+=aecalhit->ncerenkov;
	nscinttot+=aecalhit->nscintillator;
 	if(i<SCEPRINT&&ievt<SCEPRINT2) std::cout<<" hit channel is "<<aecalhit->cellID<<" in hex is "<< std::hex<< aecalhit->cellID<<std::dec<<" "<<aecalhit->energyDeposit<<" "<<aecalhit->ncerenkov<<" "<<aecalhit->nscintillator<<std::endl;




	// see ../src/DRCrystal_geo.cpp to see the assignments
	int ihitchan=aecalhit->cellID;
	int idet = (ihitchan & 0xC0)>>6;  // this assignment is made in SCEPCALConstants.xml
	int ilayer = (ihitchan & 0x38)>>3; // this is 1 for crystal and detectors, 0 for air around it
	int islice = (ihitchan & 0x07);  //   this is 1 or 4 for photodetectors, 2 for crystal
	// channels are 64 air
	//             73 77 detectors
	//            75 crystal
  //Others for the silicone gap 
	if(i<SCEPRINT&&ievt<SCEPRINT2) std::cout<<" idet,ilayer,islice is ("<<idet<<","<<ilayer<<","<<islice<<")"<<std::endl;

  // print out wavelength spectra and also timing histograms
	int ijchan=aecalhit->nbin;
 	float binsize = (aecalhit->wavelenmax-aecalhit->wavelenmin)/ijchan; //Copied from SDAction
  std::cout << "Bin size " << binsize << std::endl;
  /*int ijchan_t=aecalhit->nbin_t;
  float binsize_t = (aecalhit->timemax-aecalhit->timemin)/ijchan_t; //Copied from SDAction*/
  float kink_wavelength = 575;
  kink_cut=(kink_wavelength-aecalhit->wavelenmin)/binsize;
  //int bincut=(wavelencut-aecalhit->wavelenmin)/binsize;
  //std::cout<<bincut<<std::endl;
	//for (int j=0;j<ijchan;j++) {
	//  std::cout<<"  ncerwave["<<j<<"]="<<(aecalhit->ncerwave)[j]<<std::endl;
	//  std::cout<<"  nscintwave["<<j<<"]="<<(aecalhit->nscintwave)[j]<<std::endl;
	//}
  hchan->Fill(aecalhit->cellID);


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
    for (int j=0; j<ijchan;j++) //Wavelength bin loop
    {
     float bin_wavelength = 0.0;
     ncerwavechan[kchan].push_back((aecalhit->ncerwave)[j]);
     nscintwavechan[kchan].push_back((aecalhit->nscintwave)[j]);
     bin_wavelength= j*binsize + aecalhit->wavelenmin;
     number_of_bins[kchan].push_back(bin_wavelength); //lower wavelength of the j th bin
    }    
    
    /*for (int j_t=0; j_t<ijchan_t;j_t++) //time bin loop
    {
     ncertimechan[kchan].push_back((aecalhit->ncertime)[j_t]);
     nscinttimechan[kchan].push_back((aecalhit->nscinttime)[j_t]);
     number_of_bins_time[kchan].push_back(j_t*binsize_t + aecalhit->timemin); //lower time stamp of the j th bin
    }*/
	}
 
         // get MC truth information about individual contributions to this hit
        Contributions zxzz=aecalhit->truth;
        float hacheck=0.;
        for(size_t j=0;j<zxzz.size(); j++) {
          //      std::cout<<"testing truth truth number "<<i<<" with pdgID "<<(zxzz.at(i)).pdgID<<std::endl;
          // other member functions are trackID, deposit, time, length, x,y,z
          hacheck+=(zxzz.at(j)).deposit;
        }
        if(i<SCEPRINT2) {
	  std::cout<<"    difference between truth sum and total deposit is "<<hacheck-aecalhit->energyDeposit<<" where "<<aecalhit->energyDeposit<<" is the hit size."<<std::endl;
          if(aecalhit->energyDeposit>0) std::cout<<"      percent diff is "<<(hacheck-aecalhit->energyDeposit)/aecalhit->energyDeposit<<std::endl;
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
       //std::cout << ncertimechan[k].size() << std::endl;       
       for(int j=0;j<ncerwavechan[k].size();j++) //Right now all sizes are the same = binsize
       {
        if(ievt==0)
        {
         ncerwavechan_avg[k].push_back(ncerwavechan[k].at(j));
         nscintwavechan_avg[k].push_back(nscintwavechan[k].at(j));

         number_of_bins_avg[k].push_back(number_of_bins[k].at(j));
        }
        else if (ievt>0 &&ievt<num_evt) //i.e. the avg vector already has this size because this condition is necessarily occurring after the 'if' one above
        {
         ncerwavechan_avg[k].at(j) += ncerwavechan[k].at(j);
         nscintwavechan_avg[k].at(j) += nscintwavechan[k].at(j);        
         //Don't change the number_of_bins_avg here
        }
       }



       /*for(int j_t=0;j_t<ncertimechan[k].size();j_t++) //Right now all sizes are the same = binsize_t
       {
        if(ievt==0)
        {
         ncertimechan_avg[k].push_back(ncertimechan[k].at(j_t));
         nscinttimechan_avg[k].push_back(nscinttimechan[k].at(j_t));

         number_of_bins_time_avg[k].push_back(number_of_bins_time[k].at(j_t));
        }
        else if (ievt>0 &&ievt<num_evt) //i.e. the avg vector already has this size because this condition is necessarily occurring after the 'if' one above
        {
         ncertimechan_avg[k].at(j_t) += ncertimechan[k].at(j_t);
         nscinttimechan_avg[k].at(j_t) += nscinttimechan[k].at(j_t);        
         //Don't change the number_of_bins_time_avg here
        }
       }*/

      }
      
      std::cout << sizeof(UV_sipm_QE_x) << ' ' << sizeof(RGB_sipm_QE_x) << std::endl;
      
      /*for(int k=0;k<nchan;k++)
       {
        for(int j=0;j<ncerwavechan[k].size();j++) //Right now all sizes are the same = binsize
        {
                 
        }
       }*/
      
       
      hcEcalE->Fill(esum/1000.);
      hcEcalncer->Fill(ncertot);
      hcEcalncer0->Fill(ncerchan[0]);
      hcEcalncer1->Fill(ncerchan[1]);
      hcEcalncer2->Fill(ncerchan[2]);
      hcEcalncer3->Fill(ncerchan[3]);
      hcEcalE0->Fill(esumchan[0]);
      hcEcalE1->Fill(esumchan[1]);
      hcEcalE2->Fill(esumchan[2]);
      hcEcalE3->Fill(esumchan[3]);
      

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
      
       /*if(ievt<SCEPRINT2) std::cout<<" total energy deposit "<<esum<<std::endl;
      float check=0.;
      for( int i=0;i<nchan;i++) {
	if(ievt<SCEPRINT2) std::cout<<"esum ["<<namechan[i]<<"]="<<esumchan[i]<<std::endl;
	check+=esumchan[i];
      }
      if(ievt<SCEPRINT2) std::cout<<" check total energy desposit "<<check<<std::endl;

      if(ievt<SCEPRINT2) std::cout<<" total number of cherenkov is "<<ncertot<<std::endl;
      check=0;
      for( int i=0;i<nchan;i++) {
	if(ievt<SCEPRINT2) std::cout<<"ncerenkov ["<<namechan[i]<<"]="<<ncerchan[i]<<std::endl;
	check+=ncerchan[i];
      }
      if(ievt<SCEPRINT2) std::cout<<" check ncerenkov "<<check<<std::endl;


      if(ievt<SCEPRINT2) std::cout<<" total number of scintillator is "<<nscinttot<<std::endl;
      check=0;
      for( int i=0;i<nchan;i++) {
	if(ievt<SCEPRINT2) std::cout<<"nscintillator ["<<namechan[i]<<"]="<<nscintchan[i]<<std::endl;
	check+=nscintchan[i];
      }
      if(ievt<SCEPRINT2) std::cout<<" check nscintillator "<<check<<std::endl;*/

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
     }
     
     /*for(int j_t=0;j_t<ncertimechan_avg[i].size();j_t++)
     {
      ncertimechan_avg[i].at(j_t) = ncertimechan_avg[i].at(j_t)/num_evt;
      nscinttimechan_avg[i].at(j_t) = nscinttimechan_avg[i].at(j_t)/num_evt;
     }*/
     
     std::cout<<"esum_avg ["<<ichan[i]<<"]="<<esumchan_avg[i]<<std::endl;
     std::cout<<"nscintillator_avg ["<<ichan[i]<<"]="<<nscintchan_avg[i]<<std::endl;
     std::cout<<"ncerenkov_avg ["<<ichan[i]<<"]="<<ncerchan_avg[i]<<std::endl;
    }
    
    /*
    for(int j=0;j<ncerwavechan_avg[0].size();j++) //Right now all sizes are the same = binsize, taken the channel 0
    {
     int ent_UG11 = 91, ent_GG420 = 91;
     if (number_of_bins_avg[2].at(j) >= UV_filter_x[0] && number_of_bins_avg[2].at(j) < UV_filter_x[ent_UG11]) //within range specified for filter (which is basically the entire range over which the bins are specified)
     {
      for(int i=0;i<ent_UG11;i++)
      {
       if (number_of_bins_avg[2].at(j) >= UV_filter_x[i] && number_of_bins_avg[2].at(j) < UV_filter_x[i + 1]) //Checking what 'bin' the wavelength lies in
       filter[2].push_back((UV_filter_y[i + 1] - UV_filter_y[i]) * (number_of_bins_avg[2].at(j) - UV_filter_x[i]) / (UV_filter_x[i + 1] - UV_filter_x[i]) + UV_filter_y[i]);
      }
     }      
     else 
     {
      filter[2].push_back(0.0); //outside range
     }
     
     if (number_of_bins_avg[1].at(j) >= GG_filter_x[0] && number_of_bins_avg[1].at(j) < GG_filter_x[ent_GG420]) //within range specified for filter (which is basically the entire range over which the bins are specified)
     {
      for(int i=0;i<ent_GG420;i++)
      {
       if (number_of_bins_avg[1].at(j) >= GG_filter_x[i] && number_of_bins_avg[1].at(j) < GG_filter_x[i + 1]) //Checking what 'bin' the wavelength lies in
       filter[1].push_back((GG_filter_y[i + 1] - GG_filter_y[i]) * (number_of_bins_avg[1].at(j) - GG_filter_x[i]) / (GG_filter_x[i + 1] - GG_filter_x[i]) + GG_filter_y[i]);
      }
     }      
     else 
     {
      filter[1].push_back(0.0); //outside range
     }
     
     
     filter[0].push_back(1.0);
     filter[3].push_back(1.0);
    }
    */
     
      for(int j=0;j<ncerwavechan_avg[0].size();j++) //Right now all sizes are the same = binsize, taken the channel 0
      {
      int ent_UV = 88, ent_RGB = 88;
      if(number_of_bins_avg[1].at(j) >= UV_sipm_QE_x[0] && number_of_bins_avg[1].at(j) < UV_sipm_QE_x[ent_UV]) //i.e. if it lies within the range of sensitivity of the UV SiPM
      //1 less than actual size i.e. 23 because of usual indexing quirks (starts from 0)
      {
       for(int i=0;i<ent_UV;i++)
       {
        if (number_of_bins_avg[1].at(j) >= UV_sipm_QE_x[i] && number_of_bins_avg[1].at(j) < UV_sipm_QE_x[i + 1]) //Checking what 'bin' the wavelength lies in 
        QE[1].push_back(/*0.43 */ ((UV_sipm_QE_y[i + 1] - UV_sipm_QE_y[i]) * (number_of_bins_avg[1].at(j) - UV_sipm_QE_x[i]) / (UV_sipm_QE_x[i + 1] - UV_sipm_QE_x[i]) + UV_sipm_QE_y[i]));
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
        QE[2].push_back(/*0.325 */ ((RGB_sipm_QE_y[i + 1] - RGB_sipm_QE_y[i]) * (number_of_bins_avg[2].at(j) - RGB_sipm_QE_x[i]) / (RGB_sipm_QE_x[i + 1] - RGB_sipm_QE_x[i]) + RGB_sipm_QE_y[i]));
       }
      } 
       else 
       {
        QE[2].push_back(0.0);
       }
       
       QE[0].push_back(1.0);
       QE[3].push_back(1.0);
      }
    
    std::vector<float>wavelencut;
    for(int i=0;i<number_of_bins_avg[0].size();i++)
    {
     wavelencut.push_back(number_of_bins_avg[0].at(i)); //cuts in general are spanning over all the bins
    }

    std::cout << wavelencut.size() << std::endl;
    float ncercutchan_avg[nchan]; //Cerenkov Photons above cutoff
    float nscintcutchan_avg[nchan]; //Scintillation photons below cutoff
    std::vector<float> ncer_cut_chan_k[nchan], nscint_cut_chan_k[nchan], optimize_C_S[nchan]; //the ones that are being graphed wrt wavelength cut
     for( int i=0;i<nchan;i++)
    {
     for(int l=0;l<number_of_bins_avg[0].size();l++)
    {
     ncerwavecutchan_avg[i].push_back(0);
     nscintwavecutchan_avg[i].push_back(0);
     filter[i].push_back(0.0);
    }    
    }

    for(int k=0;k<wavelencut.size();k++)
    {
     for( int n=0;n<nchan;n++)
     {
      for(int m=0;m<number_of_bins_avg[0].size();m++)
      {
       filter[n].at(m)=0.0;
      }
     }
    
    for(int j=0;j<ncerwavechan_avg[0].size();j++) //Right now all sizes are the same = binsize, taken the channel 0
    {
     if (number_of_bins_avg[1].at(j) <= wavelencut.at(k)) //(UV SiPM side, so this is exclusively for scintillation)
     {
      filter[1].at(j)=1.0;
     }      
     else if (number_of_bins_avg[1].at(j) > wavelencut.at(k))
     {
      filter[1].at(j)=0.0;
     }
     
     if (number_of_bins_avg[2].at(j) >= wavelencut.at(k)) //(RGB SiPM side, so this is exclusively for Cerenkov)
     {
      filter[2].at(j)=1.0;
     }      
     else if (number_of_bins_avg[2].at(j) < wavelencut.at(k))
     {
      filter[2].at(j)=0.0;
     }
     
     filter[0].at(j)=1.0;
     filter[3].at(j)=1.0;
    }
    
    for( int i=0;i<nchan;i++)
    {
     for(int l=0;l<number_of_bins_avg[0].size();l++)
    {
     ncerwavecutchan_avg[i].at(l)=0.;
     nscintwavecutchan_avg[i].at(l)=0.;
    }    
     ncercutchan_avg[i] = 0.;
     nscintcutchan_avg[i] = 0.;

     for(int j=0;j<ncerwavechan_avg[i].size();j++)
     {
      /*ncerwavecutchan_avg[i].at(j) = ncerwavecutchan_avg[i].at(j)/num_evt;
      nscintwavecutchan_avg[i].at(j) = nscintwavecutchan_avg[i].at(j)/num_evt;*/
      
      ncerwavecutchan_avg[i].at(j)+=ncerwavechan_avg[i].at(j)*filter[i].at(j);
      nscintwavecutchan_avg[i].at(j)+=nscintwavechan_avg[i].at(j)*filter[i].at(j);
      
      ncerwavecutchan_avg[i].at(j) = ncerwavecutchan_avg[i].at(j)*QE[i].at(j);
      nscintwavecutchan_avg[i].at(j) = nscintwavecutchan_avg[i].at(j)*QE[i].at(j);
      
      /*ncerwavecutchan_avg[i].at(j)+=ncerwavechan_avg[i].at(j)*QE[i].at(j);
      nscintwavecutchan_avg[i].at(j)+=nscintwavechan_avg[i].at(j)*QE[i].at(j);*/
      
      ncercutchan_avg[i]+= ncerwavecutchan_avg[i].at(j);
      nscintcutchan_avg[i]+= nscintwavecutchan_avg[i].at(j);
     }
     //std::cout<<"nscintillator_avg with QE ["<<ichan[i]<<"]="<<nscintcutchan_avg[i]<<std::endl;
     //std::cout<<"ncerenkov_avg with QE ["<<ichan[i]<<"]="<<ncercutchan_avg[i]<<std::endl;     

     ncer_cut_chan_k[i].push_back(ncercutchan_avg[i]);
     nscint_cut_chan_k[i].push_back(nscintcutchan_avg[i]);
     //optimize_C_S[i].push_back(ncer_cut_chan_k[i].at(k)/nscint_cut_chan_k[i].at(k));
     optimize_C_S[i].push_back(ncer_cut_chan_k[i].at(k)/((ncerchan_avg[i])*(a/2+sqrt(nscint_cut_chan_k[i].at(k))))); //1 = a/2, here a corresponds to 2 sigma
     /*if(i==2)
     {
     std::cout << k << " " << number_of_bins_avg[i].at(k) <<" ncerenkov_avg with filter and QE ["<<ichan[i]<<"]="<<ncercutchan_avg[i] <<" nscintillator_avg with filter and QE ["<<ichan[i]<<"]="<<nscintcutchan_avg[i]<< " Optimization value "<< optimize_C_S[i].at(k)<< std::endl;
     }*/
    }
   }  
   //std::cout<<ncer_cut_chan_k[2].size()<< " " << nscint_cut_chan_k[2].size() << std::endl;
   
   /*for(int j=0;j<nscint_cut_chan_k[2].size();j++)
   {
    if(nscint_cut_chan_k[2].at(j)==0.0) //i.e. the ratio is infinity/Nan
    {
     optimize_C_S[2].at(j)=100.0; //Here 100 is a relatively large value for representation purposes
    }
   }*/
   
    float norm_cer = *max_element(ncerwavechan_avg[0].begin(), ncerwavechan_avg[0].end());
    float norm_scint = *max_element(nscintwavechan_avg[0].begin(), nscintwavechan_avg[0].end());   
    float norm_cer_det_right = *max_element(ncerwavechan_avg[2].begin(), ncerwavechan_avg[2].end());
    float norm_cer_det_left = *max_element(ncerwavechan_avg[1].begin(), ncerwavechan_avg[1].end());    
    float norm_scint_det_right = *max_element(nscintwavechan_avg[2].begin(), nscintwavechan_avg[2].end());   
    float norm_scint_det_left = *max_element(nscintwavechan_avg[1].begin(), nscintwavechan_avg[1].end());  
    float norm_S_Cut_right = *max_element(nscint_cut_chan_k[2].begin(),nscint_cut_chan_k[2].end());
    float norm_C_Cut_right = *max_element(ncer_cut_chan_k[2].begin(),ncer_cut_chan_k[2].end());
    float Large_value = *max_element(optimize_C_S[2].begin(),optimize_C_S[2].end());
    float Large_value_QE = *max_element(QE[2].begin(),QE[2].end());
    float Optimum_cut_value, Optimum_Cerenkov, Optimum_Scintillation;
    float kink_cut_value = 0.0, kink_Cerenkov = 0.0, kink_Scintillation = 0.0;
    std::cout << nscint_cut_chan_k[2].size() << std::endl;
    for(int j=0;j<nscint_cut_chan_k[2].size();j++)
    {
     if(optimize_C_S[2].at(j) == Large_value)
     {
      Optimum_cut_value = number_of_bins_avg[2].at(j);
      Optimum_Cerenkov = ncer_cut_chan_k[2].at(j);
      Optimum_Scintillation = nscint_cut_chan_k[2].at(j);
      std::cout << "Optimum cut value "<< Optimum_cut_value << " " << Optimum_Cerenkov << " " << Optimum_Scintillation << std::endl;
     }
     
     if(j == kink_cut) //Related to kink; only for PbWO4 !!!!!
     {
      kink_cut_value = number_of_bins_avg[2].at(j);
      kink_Cerenkov = ncer_cut_chan_k[2].at(j);
      kink_Scintillation = nscint_cut_chan_k[2].at(j);
      std::cout << "Kink wavelength "<< kink_cut_value << " " << kink_Cerenkov << " " << kink_Scintillation << std::endl;
     }
     
     //optimize_C_S[2].at(j) = optimize_C_S[2].at(j)/Large_value;
    }   

    for(int j=1;j<nscint_cut_chan_k[2].size()-1;j++) //For this, does not reach the end, so do not include the last and first element in the loop - quick solution
    {
     if(nscint_cut_chan_k[2].at(j) == 0 && nscint_cut_chan_k[2].at(j+1) == 0 && nscint_cut_chan_k[2].at(j-1) != 0)
     {
      std::cout << "End of scintillation spectrum " << number_of_bins_avg[2].at(j) << " " << ncer_cut_chan_k[2].at(j) << " " << nscint_cut_chan_k[2].at(j) << std::endl;
      break;
     }
    }   

    /*for(int j=0;j<nscint_cut_chan_k[2].size();j++)
   {
    if(nscint_cut_chan_k[2].at(j)==0.0) //i.e. the ratio is infinity/Nan
    {
     optimize_C_S[2].at(j)=Large_value;
    }
   }*/

    //std::cout<< "max cer is " << " " << norm_cer << endl;
    //std::cout<< "max scint is " << " " << norm_scint << endl;
    std::cout<< "max cer det right is " << " " << norm_cer_det_right << endl;
    std::cout<< "max scint det right is " << " " << norm_scint_det_right << endl;
    //std::cout<< "max cer det left is " << " " << norm_cer_det_left << endl;
    //std::cout<< "max scint det left is " << " " << norm_scint_det_left << endl;

        
    /*for( int i=0;i<nchan;i++)
    {*/
     /*for(int j=0;j<ncerwavechan_avg[0].size();j++)
     {
      ncerwavechan_avg[0].at(j) = ncerwavechan_avg[0].at(j)/norm_cer;
      nscintwavechan_avg[0].at(j) = nscintwavechan_avg[0].at(j)/norm_scint;
      ncerwavecutchan_avg[0].at(j) = ncerwavecutchan_avg[0].at(j)/norm_cer;
      nscintwavecutchan_avg[0].at(j) = nscintwavecutchan_avg[0].at(j)/norm_scint;
     }*/
    //}

     for(int j=0;j<ncerwavechan_avg[2].size();j++)
     {
      ncerwavechan_avg[2].at(j) = ncerwavechan_avg[2].at(j)*Large_value/norm_cer_det_right;
      //ncerwavecutchan_avg[2].at(j) = ncerwavecutchan_avg[2].at(j)/norm_cer_det_right/2;
      
      ncerwavechan_avg[1].at(j) = ncerwavechan_avg[1].at(j)*Large_value/norm_cer_det_left;
      //ncerwavecutchan_avg[1].at(j) = ncerwavecutchan_avg[1].at(j)/norm_cer_det_left/2;
      
      nscintwavechan_avg[1].at(j) = nscintwavechan_avg[1].at(j)*Large_value/norm_scint_det_left;
      //nscintwavecutchan_avg[1].at(j) = nscintwavecutchan_avg[1].at(j)/norm_scint_det_left/2;
      
      nscintwavechan_avg[2].at(j) = nscintwavechan_avg[2].at(j)*Large_value/norm_scint_det_right;
      //nscintwavecutchan_avg[2].at(j) = nscintwavecutchan_avg[2].at(j)/norm_scint_det_right/2;
      
      /*ncerwavechan_avg[2].at(j) = ncerwavechan_avg[2].at(j)*Large_value_QE/norm_cer_det_right;
      nscintwavechan_avg[2].at(j) = nscintwavechan_avg[2].at(j)*Large_value_QE/norm_scint_det_right;
      
      ncerwavecutchan_avg[2].at(j) = ncerwavecutchan_avg[2].at(j)*Large_value_QE/norm_cer_det_right;
      nscintwavecutchan_avg[2].at(j) = nscintwavecutchan_avg[2].at(j)*Large_value_QE/norm_scint_det_right;*/
     }
     
     for(int j=0;j<ncerwavechan_avg[0].size();j++)
     {
      nscint_cut_chan_k[2].at(j) = nscint_cut_chan_k[2].at(j)*Large_value/norm_S_Cut_right;
      ncer_cut_chan_k[2].at(j) = ncer_cut_chan_k[2].at(j)*Large_value/norm_C_Cut_right;
     }
     
    TGraph* gcerwave_gen = new TGraph(ncerwavechan_avg[0].size(),&number_of_bins_avg[0][0],&ncerwavechan_avg[0][0]);
    gcerwave_gen->SetTitle("Cerwave_gen");
    TGraph* gcerwave_det_right = new TGraph(ncerwavechan_avg[2].size(),&number_of_bins_avg[2][0],&ncerwavechan_avg[2][0]);
    gcerwave_det_right->SetTitle("Cerwave_det_right");
    TGraph* gcerwave_det_left = new TGraph(ncerwavechan_avg[1].size(),&number_of_bins_avg[1][0],&ncerwavechan_avg[1][0]);
    gcerwave_det_left->SetTitle("Cerwave_det_left");
    
     
    TGraph* gscintwave_gen = new TGraph(nscintwavechan_avg[0].size(),&number_of_bins_avg[0][0],&nscintwavechan_avg[0][0]);
    gscintwave_gen->SetTitle("Scintwave_gen");
    TGraph* gscintwave_det_left = new TGraph(nscintwavechan_avg[1].size(),&number_of_bins_avg[1][0],&nscintwavechan_avg[1][0]);
    gscintwave_det_left->SetTitle("Scintwave_det_left");
    TGraph* gscintwave_det_right = new TGraph(nscintwavechan_avg[2].size(),&number_of_bins_avg[2][0],&nscintwavechan_avg[2][0]);
    gscintwave_det_right->SetTitle("Scintwave_det_right");
    
    TGraph* gcerwave_gen_cut = new TGraph(ncerwavecutchan_avg[0].size(),&number_of_bins_avg[0][0],&ncerwavecutchan_avg[0][0]);
    gcerwave_gen_cut->SetTitle("Cerwave_gen_cut");
    TGraph* gcerwave_det_cut_right = new TGraph(ncerwavecutchan_avg[2].size(),&number_of_bins_avg[2][0],&ncerwavecutchan_avg[2][0]);
    gcerwave_det_cut_right->SetTitle("Cerwave_det_cut_right");
    TGraph* gcerwave_det_cut_left = new TGraph(ncerwavecutchan_avg[1].size(),&number_of_bins_avg[1][0],&ncerwavecutchan_avg[1][0]);
    gcerwave_det_cut_left->SetTitle("Cerwave_det_cut_left");
     
    TGraph* gscintwave_gen_cut = new TGraph(nscintwavecutchan_avg[0].size(),&number_of_bins_avg[0][0],&nscintwavecutchan_avg[0][0]);
    gscintwave_gen_cut->SetTitle("Scintwave_gen_cut");
    TGraph* gscintwave_det_cut_left = new TGraph(nscintwavecutchan_avg[1].size(),&number_of_bins_avg[1][0],&nscintwavecutchan_avg[1][0]);
    gscintwave_det_cut_left->SetTitle("Scintwave_det_cut_left");
    TGraph* gscintwave_det_cut_right = new TGraph(nscintwavecutchan_avg[2].size(),&number_of_bins_avg[2][0],&nscintwavecutchan_avg[2][0]);
    gscintwave_det_cut_right->SetTitle("Scintwave_det_cut_right");
    
    /*TGraph* gcertime_gen = new TGraph(ncertimechan_avg[0].size(),&number_of_bins_time_avg[0][0],&ncertimechan_avg[0][0]);
    gcertime_gen->SetTitle("Certime_gen");
    TGraph* gcertime_det_right = new TGraph(ncertimechan_avg[2].size(),&number_of_bins_time_avg[2][0],&ncertimechan_avg[2][0]);
    gcertime_det_right->SetTitle("Certime_det_right");*/    
    
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
    
    TGraph* filter_GG = new TGraph(filter[1].size(),&number_of_bins_avg[1][0],&filter[1][0]);
    filter_GG->SetTitle("filter for Left killMedia (GG475)");
    TGraph* filter_UV = new TGraph(filter[2].size(),&number_of_bins_avg[2][0],&filter[2][0]);
    filter_UV->SetTitle("filter for Right killMedia (UG11)");      
    
    TGraph* C_Cut = new TGraph(ncer_cut_chan_k[2].size(),&number_of_bins_avg[2][0], &ncer_cut_chan_k[2][0]);
    C_Cut->SetTitle("Distribution of integral of C above the cut");
    TGraph* S_Cut = new TGraph(nscint_cut_chan_k[2].size(),&number_of_bins_avg[2][0], &nscint_cut_chan_k[2][0]);    
    S_Cut->SetTitle("Distribution of integral of S above the cut");
    TGraph* C_S_Cut = new TGraph(optimize_C_S[2].size(),&number_of_bins_avg[2][0], &optimize_C_S[2][0]);
    C_S_Cut->SetTitle("Optimized cut value");
    
    TMultiGraph* cerscint = new TMultiGraph();
    cerscint->SetName("cer_scint");
    TCanvas *c1 = new TCanvas("c1","Graph Draw Options", 200,10,600,400); //represents coordinates of start and end points of canvas
    cerscint -> Add(gcerwave_det_right);
    //cerscint -> Add(gcerwave_det_left);    
    //cerscint -> Add(gscintwave_det_left);
    cerscint -> Add(gscintwave_det_right);
    
    //cerscint -> Add(gcerwave_det_cut_right);
    //cerscint -> Add(gcerwave_det_cut_left);    
    //cerscint -> Add(gscintwave_det_cut_left);
    //cerscint -> Add(gscintwave_det_cut_right);
    
    cerscint-> Add(C_Cut);
    cerscint-> Add(S_Cut);
    cerscint-> Add(C_S_Cut);
    
    gcerwave_det_right->SetLineColor(kRed);
    //gcerwave_det_left->SetLineColor(kRed);
    //gscintwave_det_left->SetLineColor(kOrange);
    gscintwave_det_right->SetLineColor(kOrange);
    
    //gcerwave_det_cut_right->SetLineColor(kBlack);
    //gcerwave_det_cut_left->SetLineColor(kBlack);
    //gscintwave_det_cut_left->SetLineColor(kViolet);
    //gscintwave_det_cut_right->SetLineColor(kViolet);
    
    /*cerscint -> Add(gcertime_gen);
    cerscint -> Add(gcertime_det_right);
    
    gcertime_gen->SetLineColor(kRed);
    gcertime_det_right->SetLineColor(kBlack);*/
    
    gcerwave_det_right->SetLineWidth(2);
    gscintwave_det_right->SetLineWidth(2);
    
    //gcerwave_det_cut_right->SetLineWidth(2);
    //gscintwave_det_cut_right->SetLineWidth(2);  
    C_Cut->SetLineColor(kBlack);
    C_Cut->SetLineWidth(3);
    
    S_Cut->SetLineColor(kViolet);
    S_Cut->SetLineWidth(3);
    
    C_S_Cut->SetLineColor(kGreen);
    C_S_Cut->SetLineWidth(5);
    
    //cerscint -> Add(QE_Crystal);
    //cerscint -> Add(QE_Air_out);
    //QE_Crystal->SetLineColor(kRed);
    //QE_Air_out->SetLineColor(kBlack);
    /*cerscint -> Add(QE_UV);
    QE_UV->SetLineColor(kBlue);
    QE_UV->SetLineWidth(5);*/

    /*cerscint -> Add(QE_RGB);
    QE_RGB->SetLineColor(kOrange);
    QE_RGB->SetLineWidth(5);*/
    
    /*cerscint -> Add(filter_Crystal);
    cerscint -> Add(filter_Air_out);*/
    //cerscint -> Add(filter_UV);
    //cerscint -> Add(filter_GG);
    //filter_Crystal->SetMarkerColor(kRed);
    //filter_Air_out->SetMarkerColor(kBlack);
    /*filter_UV->SetLineColor(kGreen);
    filter_UV->SetLineStyle(kDashed);    
    filter_UV->SetLineWidth(5);*/
    //filter_UV->SetMarkerColor(kBlue);
    /*filter_GG->SetLineColor(kGreen);
    filter_GG->SetLineStyle(kDashed);
    filter_GG->SetLineWidth(5);*/
    //filter_GG->SetMarkerColor(kOrange);
    
     

    cerscint -> Draw("AC"); //Draw with axes, curve
    //c1->SetLogy();
    TH1F* hist = cerscint->GetHistogram(); //Done purely for tweaking title, axes, etc.
    Double_t ymin = hist->GetMinimum();
    Double_t ymax = hist->GetMaximum();
    //The optimization cut y axis is the main one (and the Cerenkov/Scintillation distribution is the scaled one) so the ymin and ymax is for the optimization cut values

    hist->GetXaxis()->SetTitle("Wavelength (in nm)");
    hist->GetYaxis()->SetTitle("Value of optimized cut");
    hist->SetTitle("C and S produced in crystal and detected for right killMedia)");
    c1 -> BuildLegend();
    
    TGaxis *C_axis = new TGaxis(1300, ymin, 1300, ymax, ymin, ymax*norm_cer_det_right/Large_value, 510, "+L"); //Should this norm_cer_det_right*62.5 be the scaling factor for ymin as well?
    //TGaxis *C_axis = new TGaxis(1300, ymin, 1300, ymax, ymin, ymax*norm_cer_det_right/Large_value_QE, 510, "+L"); //Should this norm_cer_det_right*62.5 be the scaling factor for ymin as well?
    C_axis->Draw();
    C_axis->SetTitle("Wavelength distribution of C detected");
    C_axis->SetLineColor(kRed);
    C_axis->SetLabelColor(kRed);
    C_axis->SetTitleColor(kRed);
    C_axis->SetLabelSize(0.025);
    C_axis->SetTitleSize(0.03);
    
    TGaxis *S_axis = new TGaxis(1250, ymin, 1250, ymax, ymin, ymax*norm_scint_det_right/Large_value, 510, "+L"); //Should this norm_scint_det_right*62.5 be the scaling factor for ymin as well?
    //TGaxis *S_axis = new TGaxis(1250, ymin, 1250, ymax, ymin, ymax*norm_scint_det_right/Large_value_QE, 510, "+L"); //Should this norm_scint_det_right*62.5 be the scaling factor for ymin as well?
    S_axis->Draw();
    S_axis->SetTitle("Wavelength distribution of S detected");
    S_axis->SetLineColor(kOrange);
    S_axis->SetLabelColor(kOrange);
    S_axis->SetTitleColor(kOrange);
    S_axis->SetLabelSize(0.025);
    S_axis->SetTitleSize(0.03);
    S_axis->SetTitleOffset(-0.6);
    
    TGaxis *Cut_value = new TGaxis (Optimum_cut_value,ymin,Optimum_cut_value,ymax,1,10,0,"");
    Cut_value->SetTickLength(0);
    Cut_value->SetLineColor(kBlue);
    Cut_value->SetLineWidth(2);
    Cut_value->Draw();

    //std::string cut_photons_total = Form("%f", Optimum_Cerenkov + Optimum_Scintillation);
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << Optimum_Cerenkov + Optimum_Scintillation; //Fixed stands for fixed point precision here (no exponent), 2 is for 2 decimal places
    std::string cut_photons_total = stream.str();

    std::stringstream stream_C;
    stream_C << std::fixed << std::setprecision(2) << Optimum_Cerenkov; //Fixed stands for fixed point precision here (no exponent), 2 is for 2 decimal places
    std::string cut_photons_C = stream_C.str();
    
    std::stringstream stream_S;
    stream_S << std::fixed << std::setprecision(2) << Optimum_Scintillation; //Fixed stands for fixed point precision here (no exponent), 2 is for 2 decimal places
    std::string cut_photons_S = stream_S.str();
    
    //std::string cut_wavelength = Form("%i", Opt_value);
    std::stringstream stream1;
    stream1 << std::fixed << std::setprecision(2) << Optimum_cut_value; //Fixed stands for fixed point precision here (no exponent), 2 is for 2 decimal places
    std::string cut_wavelength = stream1.str();

    std::string photon_argument = "Total number of photons above the cut:" + cut_photons_total;
    std::string C_photon_argument = "Total number of Cerenkov photons above the cut:" + cut_photons_C;
    std::string S_photon_argument = "Total number of Scintillation photons above the cut:" + cut_photons_S;
    
    std::string wavelength_argument = "Optimum cut wavelength:" + cut_wavelength + "nm";
    std::cout << cut_photons_total << std::endl;
    std::cout << cut_wavelength << std::endl;
    
    //TPaveText *cut_photons = new TPaveText(850,0.006,1170,0.01); //For BGO
    TPaveText *cut_photons = new TPaveText(780,0.02,1170,0.03); //For PbWO4
    cut_photons->AddText(photon_argument.c_str());
    cut_photons->AddText(C_photon_argument.c_str());
    cut_photons->AddText(S_photon_argument.c_str());
    cut_photons->AddText(wavelength_argument.c_str());
    cut_photons->SetTextSize(0.02);
    cut_photons->Draw();

    //This kink drawing is only for PBWO4 !!!!!!!!!!!!!!
    
    /*TGaxis *lower_region1 = new TGaxis (550,ymin,550,ymax,1,10,0,"");
    lower_region1->SetTickLength(0);
    lower_region1->SetLineColor(kMagenta);
    lower_region1->SetLineWidth(3);
    lower_region1->Draw();
    
    TGaxis *kink_value_lower = new TGaxis (575,ymin,575,ymax,1,10,0,"");
    kink_value_lower->SetTickLength(0);
    kink_value_lower->SetLineColor(kCyan);
    kink_value_lower->SetLineWidth(3);
    kink_value_lower->Draw();
    
    TGaxis *kink_value_upper = new TGaxis (600,ymin,600,ymax,1,10,0,"");
    kink_value_upper->SetTickLength(0);
    kink_value_upper->SetLineColor(kCyan);
    kink_value_upper->SetLineWidth(3);
    kink_value_upper->Draw();*/
    
    /*int kink_value = kink_cut_value;
    std::string kink_photons_total = Form("%f", kink_Cerenkov + kink_Scintillation);
    std::string kink_wavelength = Form("%i", kink_value);
    std::string kink_photon_argument = "Total photons above the kink value:" + kink_photons_total;
    std::string kink_wavelength_argument = "Kink wavelength:" + kink_wavelength + "nm";
    std::cout << kink_photons_total << std::endl;
    
    cut_photons->AddText(kink_photon_argument.c_str());
    cut_photons->AddText(kink_wavelength_argument.c_str());
    cut_photons->SetTextSize(0.02);
    cut_photons->Draw();*/
    
        
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


