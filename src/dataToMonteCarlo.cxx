#include <TMath.h>
#include "dataToMonteCarlo.h"

bool dataToMonteCarlo(string requestedKinematic , int runNumber = 0, string targ = "", string loop = "", int alRunNumber = 0) {
  // Global ROOT settings
  gStyle->SetOptStat(0);

  //No comparison output for dummy at this time.
  if(targ == "dummy") {return false;}
  if(targ == "carbon") {return false;}
  if(targ == "chole") {return false;}
  if(targ == "coptics") {return false;}

  double beamEnergy, hsec, thetac;
  //Load in plot limits for requestedKinematic based on the kinematic setting.  x, q2, w2 for now.

  //Set kinematic limits based on central kinematics
  SetXLimits(requestedKinematic, xFocalXMin, xFocalXMax, xpFocalXMin, xpFocalXMax,
	     yFocalXMin, yFocalXMax, ypFocalXMin, ypFocalXMax, yTarXMin, yTarXMax, 
	     xpTarXMin, xpTarXMax, ypTarXMin, ypTarXMax, deltaXMin, deltaXMax,
	     w2XMin, w2XMax, q2XMin, q2XMax, xbjXMin, xbjXMax);

  /*
  cout << xFocalXMin<< " " << xFocalXMax<< " " << xpFocalXMin<< " " << xpFocalXMax<< " " <<
	     yFocalXMin<< " " << yFocalXMax<< " " << ypFocalXMin<< " " << ypFocalXMax<< " " << yTarXMin<< " " << yTarXMax<< " " << 
	     xpTarXMin<< " " << xpTarXMax<< " " << ypTarXMin<< " " << ypTarXMax<< " " << deltaXMin<< " " << deltaXMax<< " " <<
    w2XMin<< " " << w2XMax<< " " << q2XMin<< " " << q2XMax<< " " << xbjXMin<< " " << xbjXMax << endl;
  */
  //Aluminum Scale Ratio for each target
  //Apply the aluminum scale factor based on the loop in the kinematic file.
  //Needs changed to each loop.  L1 - LD2, L2 & L3 - LH2
  Float_t dummyFactor;
  static const Double_t  dFactor_loop2 = 1. / 3.789;
  static const Double_t  dFactor_loop3 = 1. / 3.789;
  static const Double_t  dFactor_loop1 = 1. / 4.063;

  bool doSubH=false, doSubD=false, doSub = false;
  if(loop == "loop1") {
    dummyFactor = dFactor_loop1;
    doSub=true;
    doSubD=true;
  } else if (loop == "loop2") {
    dummyFactor = dFactor_loop2;
    doSub=true;
    doSubH=true;
  } else if (loop == "loop3") { 
    dummyFactor = dFactor_loop3;
    doSub=true;
    doSubH=true;
  } else {
    dummyFactor = -100;
  }
  
//============================================
//   Input Kinematics Settings for MC
//============================================
/*
  string kinFile         = kinFileDir;
  fstream kinTblFile(kinFile);
  //Kin0 has a 2% offset in hsec.  
  GetTblLine(kinTblFile, kin+1);  //Changes pointer.
  string kinn, spec, name;
  kinTblFile >> kinn >> spec >> name >> beamEnergy >> hsec >> thetac;
  cout << beamEnergy << " " << hsec << " " << thetac << endl;
  setXLimits(beamEnergy, hsec, thetac);
*/
  
  // Obtain the comparison ROOT file
  cout << Form("monteCarlos/weighted/%s-%s_monte_carlo.root",requestedKinematic.c_str(),targ.c_str()) << endl;
  inputMCFile =  new TFile(Form("monteCarlos/weighted/%s-%s_monte_carlo.root",requestedKinematic.c_str(),targ.c_str()), "READ");
  mcWgtDir = (TDirectory*) (inputMCFile->FindObjectAny("mcWgtDir"));

  cout << Form("data/normalized/%s_%d.root",requestedKinematic.c_str(),runNumber) << endl;
  inputRunFile =  new TFile(Form("data/normalized/%s_%d.root",requestedKinematic.c_str(),runNumber), "READ");
  normDataDir  = (TDirectory*) (inputRunFile->FindObjectAny("dataHistogramsNorm"));
  if(!normDataDir) {
    cout << "No dataHistogramsNorm directory in TFile.  Cannot proceed.\n";
    return false;
  }

  //Set the output file
  outputComparisonFile = new TFile(Form("comps/%s/%s/%s_%s_%d.root",requestedKinematic.c_str(), targ.c_str(), requestedKinematic.c_str(), targ.c_str(),runNumber),"RECREATE");


  // Obtain the directories which contain the histograms of interest
  if(alRunNumber == 0) {
    inputAlFile =  new TFile(Form("data/normalized/%s-al.root",requestedKinematic.c_str()), "READ");
  } else {
    inputAlFile =  new TFile(Form("data/normalized/%s_%d.root",requestedKinematic.c_str(),alRunNumber), "READ");
  }
  if (!inputAlFile || inputAlFile->IsZombie()) {
    std::cerr << "No aluminum dummy to subtract.\n" << std::endl;
    doSub = false;
  }
  if(doSub) {  
    normAlDataDir  = (TDirectory*) (inputAlFile->FindObjectAny("dataHistogramsNorm"));
    if(!normAlDataDir) {
      std::cerr << "No directory in aluminum dummy file to subtract.\n" << std::endl;
      doSub = false;
    }
  }
  
  // Obtain the 1D histograms of interest
  //Monte Carlo 1D histograms
  h_xFocalMCWgt  = (TH1D*) (mcWgtDir->FindObjectAny("h_xFocalMCWgt"));
  h_xpFocalMCWgt  = (TH1D*) (mcWgtDir->FindObjectAny("h_xpFocalMCWgt"));
  h_yFocalMCWgt  = (TH1D*) (mcWgtDir->FindObjectAny("h_yFocalMCWgt"));
  h_ypFocalMCWgt  = (TH1D*) (mcWgtDir->FindObjectAny("h_ypFocalMCWgt"));
  h_yTarMCWgt  = (TH1D*) (mcWgtDir->FindObjectAny("h_yTarMCWgt"));
  h_xpTarMCWgt = (TH1D*) (mcWgtDir->FindObjectAny("h_xpTarMCWgt"));
  h_ypTarMCWgt = (TH1D*) (mcWgtDir->FindObjectAny("h_ypTarMCWgt"));
  h_deltaMCWgt = (TH1D*) (mcWgtDir->FindObjectAny("h_deltaMCWgt"));
  h_q2MCWgt    = (TH1D*) (mcWgtDir->FindObjectAny("h_q2MCWgt"));
  h_w2MCWgt    = (TH1D*) (mcWgtDir->FindObjectAny("h_w2MCWgt"));
  h_xbjMCWgt    = (TH1D*) (mcWgtDir->FindObjectAny("h_xbjMCWgt"));
  //h_wMCWgt    = (TH1D*) (mcWgtDir->FindObjectAny("h_wMCWgt"));

  //Charge and Efficiency Normalized Data histograms
  h_xFocalData = (TH1D*) (normDataDir->FindObjectAny("h_xFocalData"));
  h_xpFocalData = (TH1D*) (normDataDir->FindObjectAny("h_xpFocalData"));
  h_yFocalData = (TH1D*) (normDataDir->FindObjectAny("h_yFocalData"));
  h_ypFocalData = (TH1D*) (normDataDir->FindObjectAny("h_ypFocalData"));
  h_yTarData   = (TH1D*) (normDataDir->FindObjectAny("h_yTarData"));
  h_xpTarData  = (TH1D*) (normDataDir->FindObjectAny("h_xpTarData"));
  h_ypTarData  = (TH1D*) (normDataDir->FindObjectAny("h_ypTarData"));
  h_deltaData  = (TH1D*) (normDataDir->FindObjectAny("h_deltaData"));
  h_q2Data     = (TH1D*) (normDataDir->FindObjectAny("h_q2Data"));
  h_w2Data     = (TH1D*) (normDataDir->FindObjectAny("h_w2Data"));
  h_xbjData    = (TH1D*) (normDataDir->FindObjectAny("h_xbjData"));
  h_wData    = (TH1D*) (normDataDir->FindObjectAny("h_wData"));

  //For aluminum subtraction.
    if(doSub) {
      h_xFocalAlData = (TH1D*) (normAlDataDir->FindObjectAny("h_xFocalData"));
      h_xpFocalAlData = (TH1D*) (normAlDataDir->FindObjectAny("h_xpFocalData"));
      h_yFocalAlData = (TH1D*) (normAlDataDir->FindObjectAny("h_yFocalData"));
      h_ypFocalAlData = (TH1D*) (normAlDataDir->FindObjectAny("h_ypFocalData"));
      h_yTarAlData   = (TH1D*) (normAlDataDir->FindObjectAny("h_yTarData"));
      h_xpTarAlData  = (TH1D*) (normAlDataDir->FindObjectAny("h_xpTarData"));
      h_ypTarAlData  = (TH1D*) (normAlDataDir->FindObjectAny("h_ypTarData"));
      h_deltaAlData  = (TH1D*) (normAlDataDir->FindObjectAny("h_deltaData"));;
      h_q2AlData     = (TH1D*) (normAlDataDir->FindObjectAny("h_q2Data"));
      h_w2AlData     = (TH1D*) (normAlDataDir->FindObjectAny("h_w2Data"));
      h_xbjAlData    = (TH1D*) (normAlDataDir->FindObjectAny("h_xbjData"));
      h_wAlData    = (TH1D*) (normAlDataDir->FindObjectAny("h_wData"));    
    }

  // Obtain the 2D histos of interest
  //From Monte Carlo
  h2_xVxpFocalMCWgt  = (TH2D*) (mcWgtDir->FindObjectAny("h2_xVxpFocalMCWgt"));
  h2_xVyFocalMCWgt   = (TH2D*) (mcWgtDir->FindObjectAny("h2_xVyFocalMCWgt"));
  h2_xVypFocalMCWgt  = (TH2D*) (mcWgtDir->FindObjectAny("h2_xVypFocalMCWgt"));
  h2_xpVyFocalMCWgt  = (TH2D*) (mcWgtDir->FindObjectAny("h2_xpVyFocalMCWgt"));
  h2_xpVypFocalMCWgt = (TH2D*) (mcWgtDir->FindObjectAny("h2_xpVypFocalMCWgt"));
  h2_yVypFocalMCWgt  = (TH2D*) (mcWgtDir->FindObjectAny("h2_yVypFocalMCWgt"));

  //From Normalized Data
  h2_xVxpFocalData   = (TH2D*) (normDataDir->FindObjectAny("h2_xVxpFocalData"));
  h2_xVyFocalData    = (TH2D*) (normDataDir->FindObjectAny("h2_xVyFocalData"));
  h2_xVypFocalData   = (TH2D*) (normDataDir->FindObjectAny("h2_xVypFocalData"));
  h2_xpVyFocalData   = (TH2D*) (normDataDir->FindObjectAny("h2_xpVyFocalData"));
  h2_xpVypFocalData  = (TH2D*) (normDataDir->FindObjectAny("h2_xpVypFocalData"));
  h2_yVypFocalData   = (TH2D*) (normDataDir->FindObjectAny("h2_yVypFocalData"));

  //Currently no position checks along spectrometer in HMS
  //h2_xVyS1XMCWgt = (TH2D*) (mcWgtDir->FindObjectAny("h2_xVyS1XMCWgt"));
  //h2_xVyS1YMCWgt = (TH2D*) (mcWgtDir->FindObjectAny("h2_xVyS1YMCWgt"));
  //h2_xVyS2XMCWgt = (TH2D*) (mcWgtDir->FindObjectAny("h2_xVyS2XMCWgt"));
  //h2_xVyS2YMCWgt = (TH2D*) (mcWgtDir->FindObjectAny("h2_xVyS2YMCWgt"));
  //h2_xVyS1XData = (TH2D*) (normDataDir->FindObjectAny("h2_xVyS1XData"));
  //h2_xVyS1YData = (TH2D*) (normDataDir->FindObjectAny("h2_xVyS1YData"));
  //h2_xVyS2XData = (TH2D*) (normDataDir->FindObjectAny("h2_xVyS2XData"));
  //h2_xVyS2YData = (TH2D*) (normDataDir->FindObjectAny("h2_xVyS2YData"));
  //h2_dipoleExitMCWgt = (TH2D*) (mcWgtDir->FindObjectAny("h2_dipoleExitMCWgt"));
  //h2_dipoleExitData = (TH2D*) (normDataDir->FindObjectAny("h2_dipoleExitData"));
  //h2_caloEntrMCWgt = (TH2D*) (mcWgtDir->FindObjectAny("h2_caloEntrMCWgt"));
  //h2_caloEntrData = (TH2D*) (normDataDir->FindObjectAny("h2_caloEntrData"));
  //h2_ngcerEntrMCWgt = (TH2D*) (mcWgtDir->FindObjectAny("h2_ngcerEntrMCWgt"));
  //h2_ngcerEntrData = (TH2D*) (normDataDir->FindObjectAny("h2_ngcerEntrData"));
  //h2_hgcerEntrMCWgt = (TH2D*) (mcWgtDir->FindObjectAny("h2_hgcerEntrMCWgt"));
  //h2_hgcerEntrData = (TH2D*) (normDataDir->FindObjectAny("h2_hgcerEntrData"));


  if(doSub) {
    //Incorporate the subtraction of aluminum dummy scaled by the factor.
    //Clone the Data into the SubData histograms.
    h_xFocalSubData = (TH1D*)  h_xFocalData->Clone("h_xFocalSubData");
    h_xpFocalSubData = (TH1D*)  h_xpFocalData->Clone("h_xpFocalSubData");
    h_yFocalSubData = (TH1D*)  h_yFocalData->Clone("h_yFocalSubData");
    h_ypFocalSubData = (TH1D*)  h_ypFocalData->Clone("h_ypFocalSubData");
    h_yTarSubData = (TH1D*)  h_yTarData->Clone("h_yTarSubData");
    h_xpTarSubData = (TH1D*)  h_xpTarData->Clone("h_xpTarSubData");
    h_ypTarSubData = (TH1D*)  h_ypTarData->Clone("h_ypTarSubData");
    h_deltaSubData = (TH1D*)  h_deltaData->Clone("h_deltaSubData");
    h_q2SubData = (TH1D*)  h_q2Data->Clone("h_q2SubData");
    h_w2SubData = (TH1D*)  h_w2Data->Clone("h_w2SubData");
    h_xbjSubData = (TH1D*)  h_xbjData->Clone("h_xbjSubData");
    h_wSubData = (TH1D*)  h_wData->Clone("h_wSubData");
    //Subtract the aluminum and scale by the dummy factor.
    h_xFocalSubData->Add(h_xFocalAlData,-1.*dummyFactor);
    h_xpFocalSubData->Add(h_xpFocalAlData,-1.*dummyFactor);
    h_yFocalSubData->Add(h_yFocalAlData,-1.*dummyFactor);
    h_ypFocalSubData->Add(h_ypFocalAlData,-1.*dummyFactor);
    h_yTarSubData->Add(h_yTarAlData,-1.*dummyFactor);
    h_xpTarSubData->Add(h_xpTarAlData,-1.*dummyFactor);
    h_ypTarSubData->Add(h_ypTarAlData,-1.*dummyFactor);
    h_deltaSubData->Add(h_deltaAlData,-1.*dummyFactor);
    h_q2SubData->Add(h_q2AlData,-1.*dummyFactor);
    h_w2SubData->Add(h_w2AlData,-1.*dummyFactor);
    h_xbjSubData->Add(h_xbjAlData,-1.*dummyFactor);
    h_wSubData->Add(h_wAlData,-1.*dummyFactor);

    //Clone the AlData histograms into the scalealdata histograms
    h_xFocalScaleAlData = (TH1D*) h_xFocalAlData->Clone("h_xFocalScaleAlData");
    h_xpFocalScaleAlData = (TH1D*) h_xpFocalAlData->Clone("h_xpFocalScaleAlData");
    h_yFocalScaleAlData = (TH1D*) h_yFocalAlData->Clone("h_yFocalScaleAlData");
    h_ypFocalScaleAlData = (TH1D*) h_ypFocalAlData->Clone("h_ypFocalScaleAlData");
    h_yTarScaleAlData = (TH1D*) h_yTarAlData->Clone("h_yTarScaleAlData");
    h_xpTarScaleAlData = (TH1D*) h_xpTarAlData->Clone("h_xpTarScaleAlData");
    h_ypTarScaleAlData = (TH1D*) h_ypTarAlData->Clone("h_ypTarScaleAlData");
    h_deltaScaleAlData = (TH1D*) h_deltaAlData->Clone("h_deltaScaleAlData");
    h_q2ScaleAlData = (TH1D*) h_q2AlData->Clone("h_q2ScaleAlData");
    h_w2ScaleAlData = (TH1D*) h_w2AlData->Clone("h_w2ScaleAlData");
    h_xbjScaleAlData = (TH1D*) h_xbjAlData->Clone("h_xbjScaleAlData");
    h_wScaleAlData = (TH1D*) h_wAlData->Clone("h_wScaleAlData");
    //Scale the Aluminum for showing on plots.
    h_xFocalScaleAlData->Scale(dummyFactor);
    h_xpFocalScaleAlData->Scale(dummyFactor);
    h_yFocalScaleAlData->Scale(dummyFactor);
    h_ypFocalScaleAlData->Scale(dummyFactor);
    h_yTarScaleAlData->Scale(dummyFactor);
    h_xpTarScaleAlData->Scale(dummyFactor);
    h_ypTarScaleAlData->Scale(dummyFactor);
    h_deltaScaleAlData->Scale(dummyFactor);
    h_q2ScaleAlData->Scale(dummyFactor);
    h_w2ScaleAlData->Scale(dummyFactor);
    h_xbjScaleAlData->Scale(dummyFactor);
    h_wScaleAlData->Scale(dummyFactor);

    //In the case that aluminum subtraction is performed,
    //The histograms to plot should be the subtracted versions.
    //At this point, it would be nice to build more logic in
    //I am instead going to change the pointers for charge normalized (not subtracted yield)
    // to point to the subtracted versions.
    dataOutDir = outputComparisonFile->mkdir("normalizedData");
    dataOutDir->cd();
    h_xFocalData->Write();
    h_xpFocalData->Write();
    h_yFocalData->Write();
    h_ypFocalData->Write();
    h_yTarData->Write();
    h_xpTarData->Write();
    h_ypTarData->Write();
    h_deltaData->Write();
    h_q2Data->Write();
    h_w2Data->Write();
    h_xbjData->Write();
    h_wData->Write();
    //Now change the pointers
    h_xFocalData = h_xFocalSubData;
    h_xpFocalData = h_xpFocalSubData;
    h_yFocalData = h_yFocalSubData;
    h_ypFocalData = h_ypFocalSubData;
    h_yTarData = h_yTarSubData;
    h_xpTarData = h_xpTarSubData;
    h_ypTarData = h_ypTarSubData;
    h_deltaData = h_deltaSubData;
    h_q2Data = h_q2SubData;
    h_w2Data = h_w2SubData;
    h_xbjData = h_xbjSubData;
    h_wData = h_wSubData;
  }
   

  //Create Targ comparison canvas and fill with distributions and ratios
  c_tarComp = new TCanvas("c_tarComp", "Comparison of Target Variables", canWidth, canHeight); 
  c_tarComp->Divide(3, 2);
  c_tarComp->cd(1);
  //Fill the first canvas (6 plots).  
  //Top 3 distributions
  h_xpTarStack = new THStack("hs_xptar","Weighted MC to Data Comparison: X'_{tar}");
  drawVar(h_xpTarStack, h_xpTarMCWgt, h_xpTarData, xpTarXMin, xpTarXMax, l_xpTarComp, doSub, h_xpTarScaleAlData);

  c_tarComp->cd(2);
  h_ypTarStack = new THStack("hs_yptar","Weighted MC to Data Comparison: Y'_{tar}");
  drawVar(h_ypTarStack, h_ypTarMCWgt, h_ypTarData, ypTarXMin, ypTarXMax, l_ypTarComp, doSub, h_ypTarScaleAlData);

  c_tarComp->cd(3);
  h_deltaStack = new THStack("hs_delta","Weighted MC to Data Comparison: #delta");
  drawVar(h_deltaStack, h_deltaMCWgt, h_deltaData, deltaXMin, deltaXMax, l_deltaComp, doSub, h_deltaScaleAlData);

  //Bottom 3 data-to-mc ratios
  c_tarComp->cd(4);
  h_xpTarRatio = drawRatio(h_xpTarRatio, h_xpTarMCWgt, h_xpTarData, ln_xpTarComp, xpTarXMin, xpTarXMax);
  h_xpTarRatio->GetYaxis()->SetTitle("X'_{tar}^{data} / X'_{tar}^{MC}");
  h_xpTarRatio->SetTitle("Ratio of X'_{tar}^{data} to X'_{tar}^{MC}");

  c_tarComp->cd(5);
  h_ypTarRatio = drawRatio(h_ypTarRatio, h_ypTarMCWgt, h_ypTarData, ln_ypTarComp, ypTarXMin, ypTarXMax);
  h_ypTarRatio->GetYaxis()->SetTitle("Y'_{tar}^{data} / Y'_{tar}^{MC}");
  h_ypTarRatio->SetTitle("Ratio of Y'_{tar}^{data} to Y'_{tar}^{MC}");

  c_tarComp->cd(6);
  h_deltaRatio = drawRatio(h_deltaRatio, h_deltaMCWgt, h_deltaData, ln_deltaComp, deltaXMin, deltaXMax);
  h_deltaRatio->GetYaxis()->SetTitle("#delta_{data} / #delta_{MC}");
  h_deltaRatio->SetTitle("Ratio of #delta_{data} to #delta_{MC}");

  //Create x and y focal plots and yTar comparison canvas and fill with distributions and ratios
  c_tarComp2 = new TCanvas("c_tarComp2", "Comparison of Target Variables 2", canWidth, canHeight); 
  c_tarComp2->Divide(3, 2);
  c_tarComp2->cd(1);
  //Fill the first canvas (6 plots).  
  //Top 3 distributions
  h_yTarStack = new THStack("hs_ytar","Weighted MC to Data Comparison: Y_{tar}");
  drawVar(h_yTarStack, h_yTarMCWgt, h_yTarData, yTarXMin, yTarXMax, l_yTarComp, doSub, h_yTarScaleAlData);

  c_tarComp2->cd(2);
  h_xFocalStack = new THStack("hs_xfocal","Weighted MC to Data Comparison: X_{Foc}");
  drawVar(h_xFocalStack, h_xFocalMCWgt, h_xFocalData, xFocalXMin, xFocalXMax, l_xFocalComp, doSub, h_xFocalScaleAlData);

  c_tarComp2->cd(3);
  h_yFocalStack = new THStack("hs_yFocal","Weighted MC to Data Comparison: Y_{Foc}");
  drawVar(h_yFocalStack, h_yFocalMCWgt, h_yFocalData, yFocalXMin, yFocalXMax, l_yFocalComp, doSub, h_yFocalScaleAlData);

  //Bottom 3 data-to-mc ratios
  c_tarComp2->cd(4);
  h_yTarRatio = drawRatio(h_yTarRatio, h_yTarMCWgt, h_yTarData, ln_yTarComp, yTarXMin, yTarXMax);
  h_yTarRatio->GetYaxis()->SetTitle("Y_{tar}^{data} / Y_{tar}^{MC}");
  h_yTarRatio->SetTitle("Ratio of Y_{tar}^{data} to Y_{tar}^{MC}");

  c_tarComp2->cd(5);
  h_xFocalRatio = drawRatio(h_xFocalRatio, h_xFocalMCWgt, h_xFocalData, ln_xFocalComp, xFocalXMin, xFocalXMax);
  h_xFocalRatio->GetYaxis()->SetTitle("x_{foc}^{data} / x_{foc}^{MC}");
  h_xFocalRatio->SetTitle("Ratio of x_{foc}^{data} to x_{foc}^{MC}");

  c_tarComp2->cd(6);
  h_yFocalRatio = drawRatio(h_yFocalRatio, h_yFocalMCWgt, h_yFocalData, ln_yFocalComp, yFocalXMin, yFocalXMax);
  h_yFocalRatio->GetYaxis()->SetTitle("Y_{foc}^{data} / Y_{foc}^{MC}");
  h_yFocalRatio->SetTitle("Ratio of Y_{foc}^{data} to Y_{foc}^{MC}");


  //Create Kinematic quantity plot
  c_kinComp = new TCanvas("c_kinComp", "Comparison of Kinematic Variables", canWidth, canHeight); 
  c_kinComp->Divide(3, 2);
  c_kinComp->cd(1);
  h_xbjStack = new THStack("hs_xbj","Weighted MC to Data Comparison: X_{bj}");
  drawVar(h_xbjStack, h_xbjMCWgt, h_xbjData, xbjXMin, xbjXMax, l_xbjComp, doSub, h_xbjScaleAlData);

  c_kinComp->cd(2);
  h_q2Stack = new THStack("hs_q2","Weighted MC to Data Comparison: Q^{2}");
  drawVar(h_q2Stack, h_q2MCWgt, h_q2Data, q2XMin, q2XMax, l_q2Comp, doSub, h_q2ScaleAlData);

  c_kinComp->cd(3);
  h_w2Stack = new THStack("hs_w2","Weighted MC to Data Comparison: W^{2}");
  drawVar(h_w2Stack, h_w2MCWgt, h_w2Data, w2XMin, w2XMax, l_w2Comp, doSub, h_w2ScaleAlData);

  c_kinComp->cd(4);
  h_xbjRatio = drawRatio(h_xbjRatio, h_xbjMCWgt, h_xbjData, ln_xbjComp, xbjXMin, xbjXMax);
  h_xbjRatio->GetYaxis()->SetTitle("X_{data} / X_{MC}");
  h_xbjRatio->SetTitle("Ratio of X_{data} to X_{MC}");

  c_kinComp->cd(5);
  h_q2Ratio = drawRatio(h_q2Ratio, h_q2MCWgt, h_q2Data, ln_q2Comp, q2XMin, q2XMax);
  h_q2Ratio->GetYaxis()->SetTitle("Q^{2}_{data} / Q^{2}_{MC}");
  h_q2Ratio->SetTitle("Ratio of Q^{2}_{data} to Q^{2}_{MC}");

  c_kinComp->cd(6);
  h_w2Ratio = drawRatio(h_w2Ratio, h_w2MCWgt, h_w2Data, ln_w2Comp, w2XMin, w2XMax);
  h_w2Ratio->GetYaxis()->SetTitle("W^{2}_{data} / W^{2}_{MC}");
  h_w2Ratio->SetTitle("Ratio of W^{2}_{data} to W^{2}_{MC}");

  //Side-by-side comparisons of 2D focal plane distributions.
  c_focalComp = new TCanvas("c_focalComp", "Comparison of Focal Plane Quantites", canWidth, canHeight); 
  c_focalComp->Divide(4, 3);
  c_focalComp->cd(1); gPad->SetLogz();
  h2_xVxpFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(2); gPad->SetLogz();
  h2_xVxpFocalData->Draw("COLZ");

  c_focalComp->cd(3); gPad->SetLogz();
  h2_xVyFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(4); gPad->SetLogz();
  h2_xVyFocalData->Draw("COLZ");

  c_focalComp->cd(5); gPad->SetLogz();
  h2_xVypFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(6); gPad->SetLogz();
  h2_xVypFocalData->Draw("COLZ");

  c_focalComp->cd(7); gPad->SetLogz();
  h2_xpVyFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(8); gPad->SetLogz();
  h2_xpVyFocalData->Draw("COLZ");

  c_focalComp->cd(9); gPad->SetLogz();
  h2_xpVypFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(10); gPad->SetLogz();
  h2_xpVypFocalData->Draw("COLZ");

  c_focalComp->cd(11); gPad->SetLogz();
  h2_yVypFocalMCWgt->Draw("COLZ");

  c_focalComp->cd(12); gPad->SetLogz();
  h2_yVypFocalData->Draw("COLZ");

  //Now we will write eveything to a new output file
  //The aluminum, scaled aluminum will just be a single run for now?

  outputComparisonFile->cd();
  //Write the canvases to the main directory
  c_tarComp->SaveAs(Form("comps/%s/%s/tar/%s_%s_%d_tar.png",requestedKinematic.c_str(), targ.c_str(), requestedKinematic.c_str(), targ.c_str(), runNumber));
  c_tarComp->Write();
  c_tarComp2->SaveAs(Form("comps/%s/%s/tar2/%s_%s_%d_tar2.png",requestedKinematic.c_str(), targ.c_str(), requestedKinematic.c_str(), targ.c_str(), runNumber));
  c_tarComp2->Write();
  c_kinComp->SaveAs(Form("comps/%s/%s/kin/%s_%s_%d_kin.png",requestedKinematic.c_str(), targ.c_str(), requestedKinematic.c_str(), targ.c_str(), runNumber));
  c_kinComp->Write();
  c_focalComp->SaveAs(Form("comps/%s/%s/foc/%s_%s_%d_foc.png",requestedKinematic.c_str(), targ.c_str(), requestedKinematic.c_str(), targ.c_str(), runNumber));
  c_focalComp->Write();
  
  //Write the Histograms to their respective directories
  mcOutDir = outputComparisonFile->mkdir("mcWgt");
  mcOutDir->cd();
  h_xFocalMCWgt->Write();
  h_xpFocalMCWgt->Write();
  h_yFocalMCWgt->Write();
  h_ypFocalMCWgt->Write();
  h_yTarMCWgt->Write();
  h_xpTarMCWgt->Write();
  h_ypTarMCWgt->Write();
  h_deltaMCWgt->Write();
  h_q2MCWgt->Write();
  h_w2MCWgt->Write();
  h_xbjMCWgt->Write();
  //h_wMCWgt->Write();  

  //Add extra histograms if aluminum subtracted
  if(doSub) { 
    alumOutDir = outputComparisonFile->mkdir("normalizedAlData");
    alumOutDir->cd();
    h_xFocalAlData->Write();
    h_xpFocalAlData->Write();
    h_yFocalAlData->Write();
    h_ypFocalAlData->Write();
    h_yTarAlData->Write();
    h_xpTarAlData->Write();
    h_ypTarAlData->Write();
    h_deltaAlData->Write();
    h_q2AlData->Write();
    h_w2AlData->Write();
    h_xbjAlData->Write();
    h_wAlData->Write();
    //Subtracted histos
    subDataOutDir = outputComparisonFile->mkdir("alSubData");
    subDataOutDir->cd();
    h_xFocalSubData->Write();
    h_xpFocalSubData->Write();
    h_yFocalSubData->Write();
    h_ypFocalSubData->Write();
    h_yTarSubData->Write();
    h_xpTarSubData->Write();
    h_ypTarSubData->Write();
    h_deltaSubData->Write();
    h_q2SubData->Write();
    h_w2SubData->Write();
    h_xbjSubData->Write();
    h_wSubData->Write();
    //Write aluminum scaled histos.
    alumScaleOutDir = outputComparisonFile->mkdir("scaledNormAlData");
    alumScaleOutDir->cd();
    h_xFocalScaleAlData->Write();
    h_xpFocalScaleAlData->Write();
    h_yFocalScaleAlData->Write();
    h_ypFocalScaleAlData->Write();
    h_yTarScaleAlData->Write();
    h_xpTarScaleAlData->Write();
    h_ypTarScaleAlData->Write();
    h_deltaScaleAlData->Write();
    h_q2ScaleAlData->Write();
    h_w2ScaleAlData->Write();
    h_xbjScaleAlData->Write();
    h_wScaleAlData->Write();
  } else {
    dataOutDir = outputComparisonFile->mkdir("normalizedData");
    dataOutDir->cd();
    h_xFocalData->Write();
    h_xpFocalData->Write();
    h_yFocalData->Write();
    h_ypFocalData->Write();
    h_yTarData->Write();
    h_xpTarData->Write();
    h_ypTarData->Write();
    h_deltaData->Write();
    h_q2Data->Write();
    h_w2Data->Write();
    h_xbjData->Write();
    h_wData->Write();
  }

  inputRunFile->Close();
  if(doSub) {
    inputAlFile->Close();
  }
  inputMCFile->Close();
  outputComparisonFile->Close();

  return true;
}
