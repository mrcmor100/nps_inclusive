//#include "readutils.h"
//#include "subtraction.h"
//#include "directories.h"
//#include "radCorTbl2.h"

#include <iostream>
#include <string>
#include <TCutG.h>
#include <TF1.h>
#include <ROOT/RDataFrame.hxx>

using namespace ROOT::RDF;
using namespace std;

/*
#define D_CALO_FP 292.64
#define D_EXIT_FP -307.
#define D_S1X_FP 52.1
#define D_S1Y_FP 61.7
#define D_S2X_FP 271.4
#define D_S2Y_FP 282.4
#define D_NGCER_FP -291.700
#define D_HGCER_FP 72.600
auto calc_dipoleX = [=] (double P_tr_x, double P_tr_xp) {return P_tr_x + P_tr_xp*D_EXIT_FP; };
auto calc_dipoleY = [=] (double P_tr_y, double P_tr_yp) {return P_tr_y + P_tr_yp*D_EXIT_FP; };
auto calc_caloX = [=] (double P_tr_x, double P_tr_xp)    {return P_tr_x + P_tr_xp*D_CALO_FP; };
auto calc_caloY = [=] (double P_tr_y, double P_tr_yp)    {return P_tr_y + P_tr_yp*D_CALO_FP; };
auto calc_S1X_X = [=] (double P_dc_x, double P_dc_xp) {return P_dc_x + P_dc_xp * D_S1X_FP; };
auto calc_S1X_Y = [=] (double P_dc_y, double P_dc_yp) {return P_dc_y + P_dc_yp * D_S1X_FP; };
auto calc_S1Y_X = [=] (double P_dc_x, double P_dc_xp) {return P_dc_x + P_dc_xp * D_S1Y_FP; };
auto calc_S1Y_Y = [=] (double P_dc_y, double P_dc_yp) {return P_dc_y + P_dc_yp * D_S1Y_FP; };
auto calc_S2X_X = [=] (double P_dc_x, double P_dc_xp) {return P_dc_x + P_dc_xp * D_S2X_FP; };
auto calc_S2X_Y = [=] (double P_dc_y, double P_dc_yp) {return P_dc_y + P_dc_yp * D_S2X_FP; };
auto calc_S2Y_X = [=] (double P_dc_x, double P_dc_xp) {return P_dc_x + P_dc_xp * D_S2Y_FP; };
auto calc_S2Y_Y = [=] (double P_dc_y, double P_dc_yp) {return P_dc_y + P_dc_yp * D_S2Y_FP; };
.Define("xDipoleExit", calc_dipoleX,{"H.dc.x_fp","H.dc.xp_fp"})
.Define("yDipoleExit", calc_dipoleY,{"H.dc.y_fp","H.dc.yp_fp"})
.Define("xCaloEntr", calc_caloX,{"H.dc.x_fp","H.dc.xp_fp"})
.Define("yCaloEntr", calc_caloY,{"H.dc.y_fp","H.dc.yp_fp"})
.Define("xNGCEREntr", calc_ngcerX,{"H.dc.x_fp","H.dc.xp_fp"})
.Define("yNGCEREntr", calc_ngcerY,{"H.dc.y_fp","H.dc.yp_fp"})
.Define("xHGCEREntr", calc_hgcerX,{"H.dc.x_fp","H.dc.xp_fp"})
.Define("yHGCEREntr", calc_hgcerY,{"H.dc.y_fp","H.dc.yp_fp"})
.Define("S1X_X",calc_S1X_X, {"H.dc.x_fp","H.dc.xp_fp"})
.Define("S1X_Y",calc_S1X_Y, {"H.dc.y_fp","H.dc.yp_fp"})
.Define("S1Y_X",calc_S1Y_X, {"H.dc.x_fp","H.dc.xp_fp"})
.Define("S1Y_Y",calc_S1Y_Y, {"H.dc.y_fp","H.dc.yp_fp"})
.Define("S2X_X",calc_S2X_X, {"H.dc.x_fp","H.dc.xp_fp"})
.Define("S2X_Y",calc_S2X_Y, {"H.dc.y_fp","H.dc.yp_fp"})
.Define("S2Y_X",calc_S2Y_X, {"H.dc.x_fp","H.dc.xp_fp"})
.Define("S2Y_Y",calc_S2Y_Y, {"H.dc.y_fp","H.dc.yp_fp"})
  //Proj. to SCIN planes
  auto h2_xVyS1XData = d2.Histo2D(m_xVyS1XData, "S1X_Y", "S1X_X");
  auto h2_xVyS1YData = d2.Histo2D(m_xVyS1YData, "S1Y_Y", "S1Y_X");
  auto h2_xVyS2XData = d2.Histo2D(m_xVyS2XData, "S2X_Y", "S2X_X");
  auto h2_xVyS2YData = d2.Histo2D(m_xVyS2YData, "S2Y_Y", "S2Y_X");
  auto h2_dipoleExitData  = d2.Histo2D(m_dipoleExitData , "yDipoleExit","xDipoleExit"   );
  auto h2_caloEntrData  = d2.Histo2D(m_caloEntrData , "yCaloEntr","xCaloEntr"   );
  auto h2_ngcerEntrData  = d2.Histo2D(m_ngcerEntrData , "yNGCEREntr","xNGCEREntr"   );
  auto h2_hgcerEntrData  = d2.Histo2D(m_hgcerEntrData , "yHGCEREntr","xHGCEREntr"   );
  h2_dipoleExitData->Write();
  h2_caloEntrData->Write();
  h2_ngcerEntrData->Write();
  h2_hgcerEntrData->Write();
  h2_xVyS1XData->Write();
  h2_xVyS1YData->Write();
  h2_xVyS2XData->Write();
  h2_xVyS2YData->Write();
*/

float  mp     = 0.9382723;
float  mp2    = mp*mp;

const char* mcComparisonFile = "./dataToMC/%s_%s_%s_nominal.root";

const Int_t NBINS = 28;
Double_t edges[NBINS + 1] = {0.35,0.4,0.45,0.50,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95, 1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35, 1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1};

TDirectory *dataDir;

TH1DModel m_dummyRadScale ("h_dummyRadScale" , "Dummy Rad Scale"                                   ,  200, 0.5,1.5);
TH1DModel m_cerr       ("h_cerrEff"        ,"Cer Eff"                                               ,100,.9,1.0);      
TH1DModel m_pionC      ("h_pionC"      , "Pion Contamination; Functional Fit; Number of Entries", 200,   0, 0.2);
TH1DModel m_calEttnData ("h_calEttnData" , "Data: Calorimeter ETotTrackNorm; E/p; Number of Entries"  , 250,  0.01, 2.5);
TH1DModel m_xFocalData ("h_xFocalData" , "Data: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm"  , 100,  -40  ,  40  );
TH1DModel m_xpFocalData("h_xpFocalData", "Data: X'_{fp}; X'_{fp}; Number of Entries / 2 mrad"   , 100, -100.0, 100.0);
TH1DModel m_yFocalData ("h_yFocalData" , "Data: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm"  , 100,  -40  ,  40  );
TH1DModel m_ypFocalData("h_ypFocalData", "Data: Y'_{fp}; Y'_{fp}; Number of Entries / 2 mrad"   , 100, -100.0, 100.0);
TH1DModel m_yTarData   ("h_yTarData"   , "Data: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm", 334,  -10  ,  10  );
TH1DModel m_xpTarData  ("h_xpTarData"  , "Data: X'_{tar}; X'_{tar}; Number of Entries / 2 mrad" , 100, -100.0, 100.0);
TH1DModel m_ypTarData  ("h_ypTarData"  , "Data: Y'_{tar}; Y'_{tar}; Number of Entries / 2 mrad" , 100, -100.0, 100.0);
TH1DModel m_deltaData  ("h_deltaData"  , "Data: #delta; #delta; Number of Entries"              , 120,  -30  ,  30  );
TH1DModel m_xbjDataStatMe ("h_xbjDataStatMe"    , "Data: x_{bj}; x_{bj}; Number of Entries"             , 320,    0.0,   3.0);
TH1DModel m_q2DataStatMe ("h_q2DataStatMe"   , "Data: Q^{2}; Q^{2} (GeV^{2}); Number of Entries", 240,   0.0,  6.0);
TH1DModel m_w2DataStatMe ("h_w2DataStatMe"   , "Data: W^{2}; W^{2} (GeV^{2}); Number of Entries", 375, -10.0, 20.0);
TH1DModel m_xbjData    ("h_xbjData"    , "Data: x_{bj}; x_{bj}; Number of Entries "             , 320,    0.0,   3.0);
TH1DModel m_xbjDataMe  ("h_xbjDataMe"  , "Data: x_{bj}; x_{bj}; Number of Entries "             , 320, 0.0, 3.0);
TH1DModel m_q2Data     ("h_q2Data"     , "Data: Q^{2}; Q^{2} (GeV^{2}); Number of Entries / 0.025 GeV^{2}", 240,   0.0,  6.0);
TH1DModel m_w2Data     ("h_w2Data"     , "Data: W^{2}; W^{2} (GeV^{2}); Number of Entries / 0.050 GeV^{2}", 375, -10.0, 20.0);
TH1DModel m_q2DataMe   ("h_q2DataMe"   , "Data: Q^{2}; Q^{2} (GeV^{2}); Number of Entries / 0.025 GeV^{2}", 240,   0.0,  6.0);
TH1DModel m_w2DataMe   ("h_w2DataMe"   , "Data: W^{2}; W^{2} (GeV^{2}); Number of Entries / 0.050 GeV^{2}", 375, -10.0, 20.0);
TH1DModel m_wData      ("h_wData"      , "Data: W; W (GeV); Number of Entries / 0.050 GeV", 3000, -3.0, 3.0);
TH1DModel m_wDataMe    ("h_wDataMe"    , "Data: W; W (GeV); Number of Entries / 0.050 GeV", 3000, -3.0, 3.0);
TH1DModel m_thetaData  ("h_thetaData"  , "Data: #theta; #theta; Number of Entries / 0.01 deg"             , 200,   8.0, 18.0);
//Check on zTar in Deuterium
TH1DModel m_zTarData1  ("h_zTarData1"  , "Data: Z_{tar}; Z_{tar} (cm); Number of Entries / 1 mm", 400,  -20  ,  20  );
TH1DModel m_zTarData2  ("h_zTarData2"  , "Data: Z_{tar} X>2; Z_{tar} (cm); Number of Entries / 1 mm", 400,  -20  ,  20  );
TH1DModel m_zTarData3  ("h_zTarData3"  , "Data: Z_{tar} X>2.2; Z_{tar} (cm); Number of Entries / 1 mm", 400,  -20  ,  20  );
TH1DModel m_w2Data2    ("h_w2Data2"    , "Data: W^{2} More Bins; W^{2} (GeV^{2}); Number of Entries / 0.050 GeV^{2}", 850, -10.0, 20.0);
TH1DModel m_xbjData2   ("h_xbjData2"   , "Data: x_{bj}; x_{bj}; Number of Entries "             , 360,    0.0,   3.0);

TH2DModel m_xVyFocalData2  ("h2_xVyFocalData2"  ,"Data: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm"  ,160, -40  , 40  ,160, -40  , 40  );
TH2DModel m_xVxpFocalData ("h2_xVxpFocalData" ,"Data: X_{fp} vs. X'_{fp}; X'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm"   ,100,-100.0,100.0,160, -40  , 40  );
TH2DModel m_xVyFocalData  ("h2_xVyFocalData"  ,"Data: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm"  ,160, -40  , 40  ,160, -40  , 40  );
TH2DModel m_xVypFocalData ("h2_xVypFocalData" ,"Data: X_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm"   , 60, -60.0, 60.0,160, -40  , 40  );
TH2DModel m_xpVyFocalData ("h2_xpVyFocalData" ,"Data: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 2 mrad"   ,160, -40  , 40  ,100,-100.0,100.0);
TH2DModel m_xpVypFocalData("h2_xpVypFocalData","Data: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X'_{fp} / 2 mrad"    , 60, -60.0, 60.0,100,-100.0,100.0);
TH2DModel m_yVypFocalData ("h2_yVypFocalData" ,"Data: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; Y_{fp} (cm) / 5 mm"   , 60, -60.0, 60.0,160, -40  , 40  );
TH2DModel m_yVxpTarData   ("h2_yVxpTarData"   ,"Data: Y_{tar} vs. X'_{tar}; X'_{tar} / 2 mrad; Y_{tar} / 1 mm"    ,200,-100.0,100.0,100,  -5  ,  5  ); 
TH2DModel m_yVypTarData   ("h2_yVypTarData"   ,"Data: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; Y_{tar} / 1 mm"    ,200,-100.0,100.0,100,  -5  ,  5  ); 
TH2DModel m_xpVypTarData  ("h2_xpVypTarData"  ,"Data: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; X'_{tar} / 2 mrad", 60, -60.0, 60.0,100,-100.0,100.0);
//Sp18 specific histograms
TH2DModel m_wcDataVdelta   ("h2_wcDataVdelta"   ,"Data: W w/xpfp Cor  vs. #delta; #delta; W (GeV)"           ,428, -10,22,1000,0,3);
TH2DModel m_xpFocalDataVwc ("h2_xpFocalDataVwc" ,"Data: X'_{fp} vs. W w/xpfp Cor; X'_{fp} / 2 mrad; W (GeV)" ,480,  -60.0,  60.0,360,0.8,1.1);
TH2DModel m_ypFocalDataVwc ("h2_ypFocalDataVwc" ,"Data: Y'_{fp} vs. W w/xpfp Cor; Y'_{fp} / 2 mrad; W (GeV)" ,800, -100.0, 100.0,360,0.8,1.1);
TH2DModel m_xFocalDataVwc  ("h2_xFocalDataVwc"  ,"Data: X_{fp}  vs. W w/xpfp Cor; X_{fp}  / 2 mrad; W (GeV)" ,400,  -40  ,  40.0,360,0.8,1.1);
TH2DModel m_yFocalDataVwc  ("h2_yFocalDataVwc"  ,"Data: Y_{fp}  vs. W w/xpfp Cor; Y_{fp}  / 2 mrad; W (GeV)" ,400,  -40  ,  40.0,360,0.8,1.1);
//Strange events in Xptar vs xpfp and xfp
TH2DModel m_xpTarVxpFocalData   ("h2_xpTarVxpFocalData"   ,"Data: X'_{tar} vs. X'_{fp}; X'_{fp} / 2 mrad; X'_{tar} / 2 mrad", 100,-100.0,100.0, 200,-100.0,100.0); 
TH2DModel m_xpTarVxFocalData   ("h2_xpTarVxFocalData"   ,"Data: X'_{tar} vs. X_{fp}; X_{fp}  (cm) / 5 mm; X'_{tar} / 2 mrad", 200, -50  , 50  , 200,-100.0,100.0); 
//SCIN locatinos
TH2DModel m_xVyS1XData ("h2_xVyS1XData" ,"Data: X_{S1X} vs. Y_{S1X}; Y_{S1X} (cm) / 5 mm; X_{S1X} (cm) / 5 mm"   ,200,-100.0,100.0,200, -100  , 100  );
TH2DModel m_xVyS1YData ("h2_xVyS1YData" ,"Data: X_{S1Y} vs. Y_{S1Y}; Y_{S1Y} (cm) / 5 mm; X_{S1Y} (cm) / 5 mm"   ,200,-100.0,100.0,200, -100  , 100  );
TH2DModel m_xVyS2XData ("h2_xVyS2XData" ,"Data: X_{S2X} vs. Y_{S2X}; Y_{S2X} (cm) / 5 mm; X_{S2X} (cm) / 5 mm"   ,200,-100.0,100.0,200, -100  , 100  );
TH2DModel m_xVyS2YData ("h2_xVyS2YData" ,"Data: X_{S2Y} vs. Y_{S2Y}; Y_{S2Y} (cm) / 5 mm; X_{S2Y} (cm) / 5 mm"   ,200,-100.0,100.0,200, -100  , 100  );
TH2DModel m_dipoleExitData ("h2_dipoleExitData" ,"Tracks Projected to Dipole Exit; Y_{pos}; X_{pos}"   ,200, -100.0, 100.0,200, -100.0, 100.0);
TH2DModel m_caloEntrData ("h2_caloEntrData" ,"Tracks Projected to Calo Face; Y_{pos}; X_{pos}"          ,200, -100.0, 100.0,200, -100.0, 100.0);
TH2DModel m_ngcerEntrData ("h2_ngcerEntrData" ,"Tracks Projected to NGCER Face; Y_{pos}; X_{pos}" ,200, -100.0, 100.0,200, -100.0, 100.0);
TH2DModel m_hgcerEntrData ("h2_hgcerEntrData" ,"Tracks Projected to HGCER Face; Y_{pos}; X_{pos}" ,200, -100.0, 100.0,200, -100.0, 100.0);

TH2DModel m_wDataVdelta   ("h2_wDataVdelta"   ,"Data: W before Cor  vs. #delta; #delta; W (GeV)"             ,428,-10,22,1000,0,3);
TH2DModel m_xpFocalDataVw ("h2_xpFocalDataVw" ,"Data: X'_{fp} vs. W before Cor; X'_{fp} / 2 mrad; W (GeV)"   ,480,-60.0,60.0,360,0.8,1.1);
TH2DModel m_ypFocalDataVw ("h2_ypFocalDataVw" ,"Data: Y'_{fp} vs. W before Cor; Y'_{fp} / 2 mrad; W (GeV)"   ,800,-100.0,100.0,360,0.8,1.1);
TH2DModel m_xFocalDataVw  ("h2_xFocalDataVw"  ,"Data: X_{fp}  vs. W before Cor; X_{fp}  / 2 mrad; W (GeV)"   ,400,  -40  ,  40,360,0.8,1.1);
TH2DModel m_yFocalDataVw  ("h2_yFocalDataVw"  ,"Data: Y_{fp}  vs. W before Cor; Y_{fp}  / 2 mrad; W (GeV)"   ,400,  -40  ,  40,360,0.8,1.1);
TH2DModel m_deltaVypTarData   ("h2_deltaVypTarData"   ,"Data: #delta w/xpfp Cor  vs. Y'_{Tar}; Y'_{Tar}; #delta",100,-100,100,428, -10,22);

TH1DModel m_xbjBadX   ("h_xbjBadX"   , "Data: Poorly constructed x_{bj}; x_{bj}; Number of Entries "             , 360,    0.0,   3.0);
TH1DModel m_xbjWgtBadX("h_xbjWgtBadX"   , "Data: Poorly constructed x_{bj} with Area Weight; x_{bj}; Number of Entries "             , 360,    0.0,   3.0);
TH1DModel m_xVyBad   ("h_xVyBad"   , "Data: Bad xVy", 24,    -12.0,   12.0);
TH1DModel m_BadEvtWgt("h_xVyBad"   , "Data: Bad xVy", 100, -2, 98);
TH1DModel m_xpVypBad ("h_xVyBad"   , "Data: Bad xVy", 24,    -12.0,   12.0);


//Definitions for AddNormalization
  ifstream rptFile;
  string reportLine;
  const char* pathToNPSReports = "./data/reports/skim_NPS_HMS_report_%d_-1.report";
  //Fields to pull from report file
  TPRegexp rPreScale3Factor("\\bPs3_factor\\b");
  TPRegexp rPreScale4Factor("\\bPs4_factor\\b");
  TPRegexp rChargeBCM4A("\\bBCM4A Beam Cut Charge\\b");
  TPRegexp rChargeBCM4C("\\bBCM4C Beam Cut Charge\\b");
  TPRegexp rChargeBCM1("\\bBCM1  Beam Cut Charge\\b");
  TPRegexp rChargeBCM2("\\bBCM2  Beam Cut Charge\\b");
  TPRegexp rPs3ComputerLT("\\bPre-Scaled Ps3 HMS Computer Live Time\\b");
  TPRegexp rPs4ComputerLT("\\bPre-Scaled Ps4 HMS Computer Live Time\\b");
  TPRegexp rTrackingEff("\\bE SING FID TRACK EFFIC\\b");
  TPRegexp rCurrentBCM4A("\\bBCM4A Beam Cut Current\\b");
  TPRegexp rCurrentBCM4C("\\bBCM4C Beam Cut Current\\b");
  TPRegexp rCurrentBCM1("\\bBCM1 Beam Cut Current\\b");
  TPRegexp rCurrentBCM2("\\bBCM2 Beam Cut Current\\b");
  TPRegexp rHTrig3("\\bhTRIG3\\b");
  TPRegexp rHTrig4("\\bhTRIG4\\b");
//Not correct match of EDTM live time
  TPRegexp rPs3ElectronicLT("\\bPre-Scaled Ps3 Total Live Time (EDTM)\\b");
  TPRegexp rPs4ElectronicLT("\\bPre-Scaled Ps4 Total Live Time (EDTM)\\b");
//Patterns of int and double to match in report file
  TPRegexp rGrabDouble("\\b\\d+.\\d+\\b");
  TPRegexp rGrabInt("\\b\\d+\\b");
  TPRegexp rGrabSign("-");
//TPRegexp rGrabDouble("\\b-\\?\\d*\\.?\\d+\\b");
//TPRegexp rGrabInt("\\b-\\?\\d+\\b");

  int preScale3Factor = -10, preScale4Factor = -10;
  double bcm4aCharge = 0., bcm4cCharge = 0.,
    bcm1Charge = 0., bcm2Charge = 0.;
  double ps3ComputerLT = 0., ps4ComputerLT = 0.;
  double trackEfficiency = 0.;
  double bcm4aCurrent = 0., bcm4cCurrent = 0.,
    bcm1Current = 0., bcm2Current = 0.;
  int hTRIG3 = 0, hTRIG4 = 0;
double edtmPs3ElectronicLT = 0, edtmPs4ElectronicLT = 0;
  //Include the charge cut value
  double cutCharge = 0.;
