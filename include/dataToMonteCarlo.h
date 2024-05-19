#pragma once

// Macro to make plots comparing weighted single-arm MC events and data
// Modified for NPS experiments by Casey Morean, cmorean@jlab.org

// Input ROOT files and directores used to create output comparison file
TFile *inputMCFile; TDirectory *mcWgtDir;
TFile *inputRunFile; TDirectory *normDataDir;
TFile *inputAlFile; TDirectory *normAlDataDir;
// Output ROOT file and directories
TFile *outputComparisonFile;
TDirectory *mcOutDir, *dataOutDir, *subDataOutDir, *alumOutDir, *alumScaleOutDir;

//Monte Carlo
// Weighted 1D MC histos
TH1D *h_xFocalMCWgt, *h_xpFocalMCWgt, *h_yFocalMCWgt, *h_ypFocalMCWgt;
TH1D *h_yTarMCWgt, *h_xpTarMCWgt, *h_ypTarMCWgt, *h_deltaMCWgt;
TH1D *h_q2MCWgt, *h_w2MCWgt, *h_xbjMCWgt, *h_wMCWgt;
// Weighted 2D MC histos
TH2D *h2_xVxpFocalMCWgt, *h2_xVyFocalMCWgt, *h2_xVypFocalMCWgt;
TH2D *h2_xpVyFocalMCWgt, *h2_xpVypFocalMCWgt, *h2_yVypFocalMCWgt;

//Normalized Data
// Data 1D histos
TH1D *h_xFocalData, *h_xpFocalData, *h_yFocalData, *h_ypFocalData;
TH1D *h_yTarData, *h_xpTarData, *h_ypTarData, *h_deltaData;
TH1D *h_q2Data, *h_w2Data, *h_xbjData, *h_wData;
// Data 2D histos
TH2D *h2_xVxpFocalData, *h2_xVyFocalData, *h2_xVypFocalData;
TH2D *h2_xpVyFocalData, *h2_xpVypFocalData, *h2_yVypFocalData;

//Aluminum
//1D Al histos
TH1D *h_xFocalAlData, *h_xpFocalAlData, *h_yFocalAlData, *h_ypFocalAlData;
TH1D *h_yTarAlData, *h_xpTarAlData, *h_ypTarAlData, *h_deltaAlData;
TH1D *h_q2AlData, *h_w2AlData, *h_xbjAlData, *h_wAlData;
//2D Al Histos
TH2D *h2_xVxpFocalAlData, *h2_xVyFocalAlData, *h2_xVypFocalAlData;
TH2D *h2_xpVyFocalAlData, *h2_xpVypFocalAlData, *h2_yVypFocalAlData;

//Aluminum Subtracted, Charge Normalized Data
// Data 1D histos
TH1D *h_xFocalSubData, *h_xpFocalSubData, *h_yFocalSubData, *h_ypFocalSubData;
TH1D *h_yTarSubData, *h_xpTarSubData, *h_ypTarSubData, *h_deltaSubData;
TH1D *h_q2SubData, *h_w2SubData, *h_xbjSubData, *h_wSubData;
//No 2D histograms are aluminum subtracted.

//1D Al histos
TH1D *h_xFocalScaleAlData, *h_xpFocalScaleAlData, *h_yFocalScaleAlData, *h_ypFocalScaleAlData;
TH1D *h_yTarScaleAlData, *h_xpTarScaleAlData, *h_ypTarScaleAlData, *h_deltaScaleAlData;
TH1D *h_q2ScaleAlData, *h_w2ScaleAlData, *h_xbjScaleAlData, *h_wScaleAlData;
//No 2D histograms for the scaled down aluminum

//Data-to-MC Ratio Histograms
TH1D *h_xFocalRatio, *h_xpFocalRatio, *h_yFocalRatio, *h_ypFocalRatio;
TH1D *h_yTarRatio, *h_xpTarRatio, *h_ypTarRatio, *h_deltaRatio;
TH1D *h_q2Ratio, *h_w2Ratio, *h_xbjRatio, *h_wRatio;
//No Data-to-MC 2D Hisogram Ratios

//Position of Particals Along Spectrometer Stack
//Currently not created for the HMS
//TH2D *h2_xVyS1XData, *h2_xVyS1YData, *h2_xVyS2XData, *h2_xVyS2YData;
//TH2D *h2_xVyS1XMCWgt, *h2_xVyS1YMCWgt, *h2_xVyS2XMCWgt, *h2_xVyS2YMCWgt;
//TH2D *h2_dipoleExitData, *h2_dipoleExitMCWgt;
//TH2D *h2_caloEntrData, *h2_caloEntrMCWgt;
//TH2D *h2_ngcerEntrData, *h2_ngcerEntrMCWgt;
//TH2D *h2_hgcerEntrData, *h2_hgcerEntrMCWgt;
// Comparison canvas'
TCanvas *c_tarComp, *c_tarComp2, *c_kinComp, *c_focalComp;

// Legends
TLegend *l_xFocalComp, *l_xpFocalComp, *l_yFocalComp, *l_ypFocalComp;
TLegend *l_yTarComp, *l_xpTarComp, *l_ypTarComp, *l_deltaComp;
TLegend *l_q2Comp, *l_w2Comp, *l_xbjComp, l_wComp;
// Lines
TLine *ln_xFocalComp, *ln_xpFocalComp, *ln_yFocalComp, *ln_ypFocalComp;
TLine *ln_yTarComp, *ln_xpTarComp, *ln_ypTarComp, *ln_deltaComp;
TLine *ln_q2Comp, *ln_w2Comp, *ln_xbjComp, *ln_wComp;

//THStacks
THStack *h_xFocalStack, *h_xpFocalStack, *h_yFocalStack, *h_ypFocalStack;
THStack *h_yTarStack, *h_xpTarStack, *h_ypTarStack, *h_deltaStack;
THStack *h_q2Stack, *h_w2Stack, *h_xbjStack, *h_wStack;

// Define constants
// Canvas size paramters
Double_t canWidth   = 1280.0;
Double_t canHeight  = 720.0;
// X-axis limits, 1st line of histograms
Double_t xFocalXMin = -45.0;
Double_t xFocalXMax = 45.0;
Double_t xpFocalXMin = -45.0;
Double_t xpFocalXMax = 45.0;
Double_t yFocalXMin = -45.0;
Double_t yFocalXMax = 45.0;
Double_t ypFocalXMin = -45.0;
Double_t ypFocalXMax = 45.0;
// X-axis limits, 2ndt line of histograms
Double_t yTarXMin   = -3.5;
Double_t yTarXMax   = 3.5;
Double_t xpTarXMin  = -80.0;
Double_t xpTarXMax  = 80.0;
Double_t ypTarXMin  = -40.0;
Double_t ypTarXMax  = 40.0;
Double_t deltaXMin  = -9.0;
Double_t deltaXMax  = 9.0;
// X-axis limits, 3rd line of histograms
//Ideally the setXLimits method sets kinematic limits
Double_t q2XMin     = 2.;
Double_t q2XMax     = 6.;
Double_t w2XMin     = 0.0;
Double_t w2XMax     = 16;
Double_t xbjXMin    = 0.15;
Double_t xbjXMax    = 0.7;
Double_t wXMin    = -4.0;
Double_t wXMax    = 8;

// Ratio limits
static const Double_t idealRatio = 1.0;
static const Double_t ratioMin   = 0.5;
static const Double_t ratioMax   = 1.5;

string limitsFilename = "kinematics/kinSettingsAllLimits.dat";

bool SetXLimits(const std::string& searchString, Double_t &var0, Double_t &var1, Double_t &var2, Double_t &var3, 
		Double_t &var4, Double_t &var5, Double_t &var6, Double_t &var7, Double_t &var8, Double_t &var9, 
		Double_t &var10, Double_t &var11, Double_t &var12, Double_t &var13, Double_t &var14, Double_t &var15, 
		Double_t &var16, Double_t &var17, Double_t &var18, Double_t &var19, Double_t &var20, Double_t &var21) {
  std::fstream file(limitsFilename, std::ios::in);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file." << std::endl;
    return false; // Exit with error code if file cannot be opened
  }
  
  std::string line;
  // Move back to the beginning of the file
  file.seekg(0, std::ios::beg);
  
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string token;
    // Get the first token which should be the string we're searching for
    if (iss >> token) {
      if (token == searchString) {
	// If the string matches, extract the other variables from the line
	if (!(iss >> var0 >> var1 >> var2 >> var3 >> 
	      var4 >> var5 >> var6 >> var7 >> var8 >> var9 >> 
	      var10 >> var11 >> var12 >> var13 >> var14 >> var15 >> 
	      var16 >> var17 >> var18 >> var19 >> var20 >> var21)) {
	  // Handle error if the line doesn't have enough data
	  std::cerr << "Error: Line format incorrect." << std::endl;
	  return false;
	}
	return true; // Found the string and extracted the data
      }
    }
  }
  // String was not found in the file
  std::cerr << "Error: String not found." << std::endl;
  return false;
}

void setXLimits(double E, double Ep, double th) {
  //Not right.
  double d2r = 3.14159 / 180.;
  double thetaCentralRadians = th * d2r;
  double mp  = 0.9383;
  double thXMin, thXMax;
  thXMin = TMath::ACos(TMath::Cos(thetaCentralRadians-ypTarXMax)*TMath::Cos(xpTarXMin));
  thXMax = TMath::ACos(TMath::Cos(thetaCentralRadians-ypTarXMin)*TMath::Cos(xpTarXMax));
  q2XMin     = 4*E*(Ep+(deltaXMin/100.))*pow(sin(thXMin*d2r/2),2);
  q2XMax     = 4*E*(Ep+(deltaXMax/100.))*pow(sin(thXMax*d2r/2),2);
  w2XMin     = mp*mp + 2*mp*(E-Ep+(deltaXMin/100.)) - q2XMax - 1;
  w2XMax     = mp*mp + 2*mp*(E-Ep+(deltaXMax/100.)) - q2XMin + 2;
  xbjXMin    = q2XMin / (2*mp*(E-Ep*(deltaXMin/100.)));
  xbjXMax    = q2XMax / (2*mp*(E-Ep*(deltaXMax/100.)));
  return;
}

void drawVar(THStack* hs, TH1D* mc, TH1D* data, double xMin, double xMax, TLegend* leg, bool doSubCryo=false, TH1D* al = NULL) {
  int mc_lc = 8; int mc_fc = 8; int mc_fs = 3008;
  int data_lc = 4; int data_ms = 22;
  int al_mc = kMagenta+1; int al_ms = 23;
  mc->SetLineColor(mc_lc);
  mc->SetFillStyle(mc_fs);
  mc->SetFillColor(mc_fc);
  hs->Add(mc,"hist");
  data->SetMarkerColor(data_lc);
  data->SetMarkerStyle(data_ms);
  hs->Add(data);
  if(doSubCryo) {
    al->SetMarkerColor(al_mc);
    al->SetMarkerStyle(al_ms);
    hs->Add(al);
  }
  leg = new TLegend(0.70, 0.70, 0.9, 0.9);
  leg->AddEntry(mc, "Weighted MC", "F");
  leg->AddEntry(data, "Data", "EP");
  if(doSubCryo) {
    leg->AddEntry(al, "Aluminum", "EP");
  }
  hs->Draw("nostack");
  hs->GetXaxis()->SetRangeUser(xMin, xMax);

  leg->Draw();
  return;
}

TH1D* drawRatio(TH1D* ratio, TH1D* mc, TH1D* data, TLine* ln, double xMin, double xMax) {
  const char* name = data->GetName();
  int r_mstyle = 22; float r_msize = 1.25;
  int ln_style = 9; int ln_width = 2; int ln_color=38;
  ratio = dynamic_cast <TH1D*> (data->Clone(name));
  ratio->Divide(mc);
  ratio->SetMarkerStyle(r_mstyle);
  ratio->SetMarkerSize(r_msize);
  ratio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  ratio->GetXaxis()->SetRangeUser(xMin, xMax);

  ratio->Draw();

  ln = new TLine(xpTarXMin, idealRatio, xpTarXMax, idealRatio);
  ln->SetLineStyle(ln_style);
  ln->SetLineWidth(ln_width);
  ln->SetLineColor(ln_color);
  ln->Draw();
  return ratio;
}
