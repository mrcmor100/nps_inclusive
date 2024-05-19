#include <iostream>
#include <string>

#include "TGraph2D.h"
#include <TCutG.h>
#include <TF1.h>
#include <ROOT/RDataFrame.hxx>

#include "weightMC.h"
#include "processRun.h"
#include "extractData.h"

using namespace ROOT::RDF;

int weightManyMC(string requestedKinematic) {
  //Now do all the targets for this kinematic
  std::set<std::string> targs;
  if(!GetTargetsInKinematic(requestedKinematic, targs)) {
    std::cerr << "Issues with targets in this kinematic." << requestedKinematic << std::endl;
      return -5; // Or handle error as needed
  }
  for(string targ : targs) {
    if(targ == "ld2" || targ == "lh2") {
      cout << "Currently no target " << targ << " in MC.\n";
      weightMC(requestedKinematic, targ);
    } else {
      cout << "Currently no target " << targ << " in MC.\n";
    }
  }
  return 0;
}

int weightMC(string requestedKinematic, string target) {

//============================================
//   Input Kinematics Settings for MC
//============================================
  std::fstream kinFile(kinematicFilename, std::ios::in);
  if (!kinFile.is_open()) {
    std::cerr << "Error: Could not open file." << std::endl;
    return -1; // Exit with error code if file cannot be opened
  }

// Search for "requestedKinematic" and extract data into spectrometer, eBeam, eCentral, thetaCentral
  if (!GetKinematicInfoByString(kinFile, requestedKinematic, spectrometer, eBeam, eCentral, thetaCentral)) {
    kinFile.close();
    return -2; // Exit with a custom error code if string not found or other error

  }
  cout << "=============================================\n";
  cout << "           Beam and Kineamtic Info           \n";
  cout << "=============================================\n";
  cout << " Beam Energy: " << eBeam << endl;
  cout << " Central Momentum: " << eCentral << endl;
  cout << " Central Angle (deg): " << thetaCentral << endl;
  cout << endl << endl;
//Convert the central angel to radians
  double thetaCentralRadians = thetaCentral * TMath::DegToRad();

//============================================
//      Input file directories and names
//============================================
  const char* mcRawInputForm = "./monteCarlos/infiles/%s-%s.inp";
  string mcRawInputFilename = Form(mcRawInputForm, requestedKinematic.c_str(), target.c_str());
  const char* mcRawForm = "./monteCarlos/worksim/%s-%s.root";
  string mcRawFiles      = Form(mcRawForm, requestedKinematic.c_str(),target.c_str());
  cout << "=============================================\n";
  cout << "           Monte Carlo Files           \n";
  cout << "=============================================\n";
  cout << "Monte Carlo Input File: " << mcRawInputFilename << endl;
  cout << "Monte Carlo Raw Data File(s): " << mcRawFiles << endl;
  cout << endl << endl;

//============================================
//     input to raw monte carlo
//============================================
  fstream mcInputFile(mcRawInputFilename, std::ios::in);
  float deltaDown = GetInpVar(mcInputFile, 8);
  float deltaUp   = GetInpVar(mcInputFile, 9);
  float xpDown = GetInpVar(mcInputFile, 12);
  float xpUp   = GetInpVar(mcInputFile, 13);
  float ypDown = GetInpVar(mcInputFile, 10);
  float ypUp   = GetInpVar(mcInputFile, 11);

  cout << "=============================================\n";
  cout << "            Phase Space Limits in MC         \n";
  cout << "=============================================\n";
  cout << " Delta: deldown = " << deltaDown << endl;
  cout << " Delta: delup   = " << deltaUp << endl;
  cout << " xpTar: dxpdown = " << xpDown << endl;
  cout << " xpTar: dxpup   = " << xpUp << endl;
  cout << " ypTar: dypdown = " << ypDown << endl;
  cout << " ypTar: dypup   = " << ypUp << endl;
  cout << endl << endl;

  double deltaEp = (deltaUp - deltaDown) / 100. * eCentral;
  double deltaXp = (xpUp - xpDown);
  double deltaYp = (ypUp - ypDown);
  double phaseSpace = deltaXp * deltaYp * deltaEp / 1000.0;
  cout << "=============================================\n";
  cout << "              Phase Space                    \n";
  cout << "=============================================\n";
  cout << "dE': " << deltaEp << endl;
  cout << "dx': " << deltaXp << endl;
  cout << "dy': " << deltaYp << endl;
  cout << "Phase Space: " <<  phaseSpace << endl;
  cout << endl << endl;

//============================================
//       Read in Target data from file
//============================================
  std::fstream targetInfoFile(targetInfoFilename, std::ios::in);
  if (!targetInfoFile.is_open()) {
    std::cerr << "Error: Could not open file: " << std::endl;
    return -1; // Exit with error code if file cannot be opened
  }

  int targetID, atomicNumber;
  float density, targetLength, atomicMass, aerialDensity, dataLuminosity, cryoShrinkFactor;
// Search for 'target' and extract data into spectrometer, eBeam, eCentral, thetaCentral
  if (!GetTargetInfoByString(targetInfoFile, target, targetID, density,
			     targetLength, atomicMass, atomicNumber, aerialDensity,
			     dataLuminosity, cryoShrinkFactor)) {
    targetInfoFile.close();
    return -2; // Exit with a custom error code if string not found or other error
  }

  dataLuminosity = dataLuminosity * cryoShrinkFactor;

  cout << "=============================================\n";
  cout << "              Target Information             \n";
  cout << "=============================================\n";
  cout << " Target:    " << target << endl;
  cout << " Target ID: " << targetID << endl;
  cout << " Target Density:  " << density << endl;
  cout << " Aerial Density:  " << aerialDensity << endl;
  cout << " Target Length:   " << targetLength << endl;
  cout << " Atomic Mass:     " << atomicMass << endl;
  cout << " Atomic Number:   " << atomicNumber << endl;
  if(std::strcmp(target.c_str(),"lh2")==0) {
    cout << " Cryo Correction: " << cryoShrinkFactor << endl;}
  else if(std::strcmp(target.c_str(),"ld2")==0) {
    cout << " Cryo Correction: " << cryoShrinkFactor << endl;}
  cout << " Data Luminosity: " << dataLuminosity << endl;
  cout << endl << endl;


//==============================================
//  Interpolation from a 2D graph using root.
//     and put them in a TGraph2D().
//==============================================
 TGraph2D* xt = new TGraph2D();
 TGraph2D* rt = new TGraph2D();
 TGraph2D* cst = new TGraph2D();

//variables to hold elements in each line
 double t_e, t_ep, t_theta, t_xbj, t_Q2, t_w2;
 double t_amu, t_radCorrFactor, t_ccor;
 double Sig_Born, Sig_Born_In, Sig_Born_QE,
   Sig_Rad, Sig_Rad_EL, Sig_Rad_QE, Sig_Rad_DIS;

 int i = 0;
 string line;
 //Three tables, one for each beam energy:
 string csTableFile = Form("crossSectionModel/%s-%s.out",requestedKinematic.c_str(),target.c_str());

 ifstream crossFile(csTableFile, std::ios::in);
//Open file

 cout << "Reading from file: " << csTableFile << endl;
  int nSkip = 2000;
  cout << "=============================================\n";
  cout << "            Cross Section Table              \n";
  cout << "=============================================\n";
  cout << "Every entry that is a multiple of " << nSkip << ":\n";
  cout << " E'\tTheta\tX\tSigRad\n";
  
  bool setAnyPts = false;
  //Skip first line
  while(getline(crossFile, line)) {
    if(i == 0) {
      i+=1; continue;
    }
    istringstream ss(line);
    i+=1;

    //Read from cross section table (output from rc-externals)
    ss >> t_e >> t_ep >> t_theta  >> t_xbj >> t_Q2 >>
      Sig_Born >> Sig_Born_In >> Sig_Born_QE >> Sig_Rad >>
      Sig_Rad_EL >> Sig_Rad_QE >> Sig_Rad_DIS >> t_ccor;

    //Calculate t_w2, not in table originally.
    t_w2 = mp2 + 2.*mp*(t_e-t_ep) - t_Q2;

    if(i%nSkip == 0) {
      cout << t_w2 << "\t" << t_theta << "\t" << t_xbj << "\t" << Sig_Rad  << endl;
    }
    
    setAnyPts = true;
    cst->SetPoint(i, t_w2, t_theta, Sig_Rad);
  }
  if(!setAnyPts) {
    cout << "Error! No cross-section!\n\n";
    return -1;
  }

  auto calc_dipoleX = [=] (float P_tr_x, float P_tr_xp) {return P_tr_x + P_tr_xp*D_EXIT_FP; };
  auto calc_dipoleY = [=] (float P_tr_y, float P_tr_yp) {return P_tr_y + P_tr_yp*D_EXIT_FP; };
  auto calc_ngcerX = [=] (float P_tr_x, float P_tr_xp) {return P_tr_x + P_tr_xp*D_NGCER_FP; };
  auto calc_ngcerY = [=] (float P_tr_y, float P_tr_yp) {return P_tr_y + P_tr_yp*D_NGCER_FP; };
  auto calc_hgcerX = [=] (float P_tr_x, float P_tr_xp) {return P_tr_x + P_tr_xp*D_HGCER_FP; };
  auto calc_hgcerY = [=] (float P_tr_y, float P_tr_yp) {return P_tr_y + P_tr_yp*D_HGCER_FP; };
  auto calc_caloX = [=] (float P_tr_x, float P_tr_xp)    {return P_tr_x + P_tr_xp*D_CALO_FP; };
  auto calc_caloY = [=] (float P_tr_y, float P_tr_yp)    {return P_tr_y + P_tr_yp*D_CALO_FP; };
  auto calc_S1X_X = [=] (float P_dc_x, float P_dc_xp) {return P_dc_x + P_dc_xp * D_S1X_FP;};
  auto calc_S1X_Y = [=] (float P_dc_y, float P_dc_yp) {return P_dc_y + P_dc_yp * D_S1X_FP;};
  auto calc_S1Y_X = [=] (float P_dc_x, float P_dc_xp) {return P_dc_x + P_dc_xp * D_S1Y_FP;};
  auto calc_S1Y_Y = [=] (float P_dc_y, float P_dc_yp) {return P_dc_y + P_dc_yp * D_S1Y_FP;};
  auto calc_S2X_X = [=] (float P_dc_x, float P_dc_xp) {return P_dc_x + P_dc_xp * D_S2X_FP;};
  auto calc_S2X_Y = [=] (float P_dc_y, float P_dc_yp) {return P_dc_y + P_dc_yp * D_S2X_FP;};
  auto calc_S2Y_X = [=] (float P_dc_x, float P_dc_xp) {return P_dc_x + P_dc_xp * D_S2Y_FP;};
  auto calc_S2Y_Y = [=] (float P_dc_y, float P_dc_yp) {return P_dc_y + P_dc_yp * D_S2Y_FP;};

  auto rad2mrad   = [](float pFocalData) {return pFocalData*1000.0;};
  auto hseCalc = [=] (float deltai) {return eCentral * (1.0 + deltai / 100.);};
  auto hseCalc2 = [=] (double deltai) {return eCentral*(1.0 + deltai / 100.);};
  auto thetainiCalc = [=] (float yptari, float xptari) {
    //Check if SHMS / HMS
    return TMath::ACos(TMath::Cos(thetaCentralRadians-yptari)*TMath::Cos(xptari));
  };
  auto hsthetaCalc = [=] (float yptar, float xptar)   {
    //Check if SHMS / HMS
    return TMath::ACos(TMath::Cos(thetaCentralRadians-yptar)*TMath::Cos(xptar));
  };
  auto sin2Calc = [] (double thetaini) {return TMath::Power(TMath::Sin(thetaini/2.), 2);};
  auto nuCalc = [=] (double eSpectrometer) {return eBeam - eSpectrometer;};
  auto q2Calc = [=] (double eSpectrometer, double sin2) {return 4.0 * eSpectrometer * eBeam * sin2;};
  auto w2Calc = [=] (double nu, double Q2) {return mp2 + 2.*mp*nu - Q2;};
  auto xbjCalc = [=] (double Q2, double nu) {return Q2 / (2*mp*nu);};
  auto jacobianCalc = [] (float psxptari, float psyptari) {return 1./TMath::Power((1+TMath::Power(psxptari,2)+TMath::Power(psyptari,2)),(3/2));};
  auto csbCalc = [=] (double thetaini, double hsev) {
    return 0.0; //Do not use now.
    double csb_cx = 0;
    Double_t p0=-2.09 * thetaini*180./TMath::Pi() +12.47;
    Double_t p1=0.2 * thetaini*180./TMath::Pi() -0.6338;
    csb_cx=exp(p0)*(exp(p1*(eBeam-hsev))-1.);
    if(std::strcmp(target.c_str(),"lh2")==0) return csb_cx;
    if(std::strcmp(target.c_str(),"ld2")==0) return csb_cx=2*csb_cx;
    if(std::strcmp(target.c_str(),"c12")==0) return csb_cx=12*csb_cx*19.32/((890.4+769.1)/2); //need to add. (use rad length)
    return 0.0;
  };
  
  auto deltaCor = [=] (float dp, float xpfp, float ypfp) {
    return (double) dp;
  };

  //Open the unweighted root tree files. chain with 1411
  ROOT::EnableImplicitMT(8);
  //cout << ROOT::GetImplicitMTPoolSize() << endl;
  ROOT::RDataFrame d("h1", mcRawFiles); // Interface to TTree and TChain with ALL monte-carlo files '*'

  auto d2 = d.Define("delta" , deltaCor, {"hsdelta","hsxpfp","hsypfp"})
    .Filter("stop_id == 0")
    .Filter("abs(hsxptar) < 0.100")
    .Filter("abs(hsyptar) < 0.050")
    .Filter("abs(hsytar) <= 10")
    .Filter("hsdelta < 8. && hsdelta > -8.");

  //Add releveant branches
  d2 = d2.Define("hsev", hseCalc, {"hsdeltai"})
         .Define("hse", hseCalc2, {"delta"})
         .Define("xptari_mrad", rad2mrad, {"hsxptari"})
         .Define("yptari_mrad", rad2mrad, {"hsyptari"})
         .Define("xptar_mrad", rad2mrad, {"hsxptar"})
         .Define("yptar_mrad", rad2mrad, {"hsyptar"})
         .Define("xpfp_mrad", rad2mrad, {"hsxpfp"})
         .Define("ypfp_mrad", rad2mrad, {"hsypfp"})
         .Define("thetaini", thetainiCalc, {"hsyptari","hsxptari"})
         .Define("hstheta", hsthetaCalc, {"hsyptar","hsxptar"})
         .Define("sin2_i", sin2Calc , {"thetaini"})
         .Define("sin2_r", sin2Calc , {"hstheta"})
         .Define("nu_i", nuCalc , {"hsev"})
         .Define("nu_r", nuCalc , {"hse"})
         .Define("q2_i", q2Calc, {"hsev","sin2_i"})
         .Define("q2_r", q2Calc, {"hse","sin2_r"})
         .Define("w2_i", w2Calc, {"nu_i","q2_i"})
         .Define("w2_r", w2Calc, {"nu_r","q2_r"})
         .Define("xbj_i", xbjCalc, {"q2_i","nu_i"})
         .Define("xbj_r", xbjCalc, {"q2_r","nu_r"})
         .Define("csb_cx",csbCalc,{"thetaini","hsev"})
         .Define("phasespcor", jacobianCalc, {"hsxptari","hsyptari"})
         .Define("xDipoleExit", calc_dipoleX,{"hsxfp","hsxpfp"})
         .Define("yDipoleExit", calc_dipoleY,{"hsyfp","hsypfp"})
         .Define("xCaloEntr", calc_caloX,{"hsxfp","hsxpfp"})
         .Define("yCaloEntr", calc_caloY,{"hsyfp","hsypfp"})
         .Define("xNGCEREntr", calc_ngcerX,{"hsxfp","hsxpfp"})
         .Define("yNGCEREntr", calc_ngcerY,{"hsyfp","hsypfp"})
         .Define("xHGCEREntr", calc_hgcerX,{"hsxfp","hsxpfp"})
         .Define("yHGCEREntr", calc_hgcerY,{"hsyfp","hsypfp"})
         .Define("S1X_X",calc_S1X_X, {"hsxfp","hsxpfp"})
         .Define("S1X_Y",calc_S1X_Y, {"hsyfp","hsypfp"})
         .Define("S1Y_X",calc_S1Y_X, {"hsxfp","hsxpfp"})
         .Define("S1Y_Y",calc_S1Y_Y, {"hsyfp","hsypfp"})
         .Define("S2X_X",calc_S2X_X, {"hsxfp","hsxpfp"})
         .Define("S2X_Y",calc_S2X_Y, {"hsyfp","hsypfp"})
         .Define("S2Y_X",calc_S2Y_X, {"hsxfp","hsxpfp"})
         .Define("S2Y_Y",calc_S2Y_Y, {"hsyfp","hsypfp"});

  if(std::strcmp(target.c_str(),"lh2")==0) {
    d2 = d2.Filter("w2_r > 0.4");
  }

  auto bornCalc = [&] (double hse, double thetaini) {return cst->Interpolate(hse, TMath::RadToDeg()*thetaini);};

  auto bornEMCCalc = [&] (double w2, double thetaini) {
    double cs = (cst->Interpolate(w2, TMath::RadToDeg()*thetaini));
    //double cs = 0.;
    //cout << "cs failed\n";
    return cs != 0 ? cs :  cs;
    //return cs;
  };

  auto ngenCutBorn = [=] (double born) {
    return born > 0.0;
  };
  auto ngenCut = [=] (float yptari, float xptari, float deltai) {
    return (xptari < xpUp/1000. && xptari > xpDown/1000. && yptari < ypUp/1000. && 
	    yptari > ypDown / 1000. && deltai > deltaDown && deltai < deltaUp);
  };

  auto d3 = d2.Define("born", bornEMCCalc, {"w2_i","thetaini"});

  cout << "=============================================\n";
  cout << "            Monte Carlo Events               \n";
  cout << "=============================================\n";
  auto ngenPassed = d3.Filter(ngenCut,{"hsyptari","hsxptari","hsdeltai"}).Count();
  cout << *ngenPassed << " events made it to the detector\n  within the phase space.\n";
  auto ngenPassedBorn = d3.Filter(ngenCutBorn,{"born"}).Count();
  cout << *ngenPassedBorn << " events passed the\n  born cut and are in phase space.\n";
  auto ngen = d.Filter(ngenCut,{"hsyptari","hsxptari","hsdeltai"}).Count();
  cout << *ngen << " events were generated in the\n  gen limits of  the monte-carlo.\n";
  cout << endl << endl;

  double fract;
  fract = dataLuminosity * phaseSpace / *ngen / 1000.00;
  cout << "=============================================\n";
  cout << "           Monte Carlo Weighting             \n";
  cout << "=============================================\n";
  cout << " Fract: " << fract << endl;
  cout << endl << endl;

  auto mcWgtCalc = [=] (double born, double phasespcor, double csb_cs, double xbj) {
    if(born > 0. && xbj >= 2.2) {
      return 0.0;
    }
    return born==0 ? 0. : (born+csb_cs)*phasespcor*fract;
  };
  const char* mcComparisonFile = "monteCarlos/weighted/%s-%s_monte_carlo.root";
  cout << "=============================================\n";
  cout << "                   Output                    \n";
  cout << "=============================================\n";
  cout << " Writing root file to: \n";
  cout << Form(mcComparisonFile, requestedKinematic.c_str(), target.c_str()) << endl;
  cout << endl << endl;

  d3 = d3.Define("mcWgt", mcWgtCalc, {"born","phasespcor","csb_cx","xbj_r"});
  
  TFile *compFile;
  compFile  = new TFile(Form(mcComparisonFile, requestedKinematic.c_str(), target.c_str()), "UPDATE");
  //Define all the histograms
  auto h_xFocalMCWgt  = d3.Histo1D(m_xFocalMCWgt,"hsxfp","mcWgt");
  auto h_xpFocalMCWgt = d3.Histo1D(m_xpFocalMCWgt,"xpfp_mrad","mcWgt");
  auto h_yFocalMCWgt  = d3.Histo1D(m_yFocalMCWgt,"hsyfp","mcWgt");
  auto h_ypFocalMCWgt = d3.Histo1D(m_ypFocalMCWgt,"ypfp_mrad","mcWgt");
  auto h_yTarMCWgt    = d3.Histo1D(m_yTarMCWgt,"hsytar","mcWgt");      //Reconstructed
  auto h_xpTarMCWgt   = d3.Histo1D(m_xpTarMCWgt,"xptar_mrad","mcWgt"); //Reconstructed
  auto h_ypTarMCWgt   = d3.Histo1D(m_ypTarMCWgt,"yptar_mrad","mcWgt"); //Reconstructed
  auto h_deltaMCWgt   = d3.Histo1D(m_deltaMCWgt,"hsdelta","mcWgt");    //Reconstructed
  auto h_q2MCWgt      = d3.Histo1D(m_q2MCWgt,"q2_r","mcWgt");
  auto h_w2MCWgt      = d3.Histo1D(m_w2MCWgt,"w2_r","mcWgt");
  auto h_xbjMCWgt     = d3.Histo1D(m_xbjMCWgt,"xbj_r","mcWgt");
  // Weighted 2D MC histos
  auto h2_xVxpFocalMCWgt = d3.Histo2D(m_xVxpFocalMCWgt,"xpfp_mrad","hsxfp","mcWgt");
  auto h2_xVyFocalMCWgt  = d3.Histo2D(m_xVyFocalMCWgt,"hsyfp","hsxfp","mcWgt");
  auto h2_xVypFocalMCWgt = d3.Histo2D(m_xVypFocalMCWgt,"ypfp_mrad","hsxfp","mcWgt");
  auto h2_xpVyFocalMCWgt = d3.Histo2D(m_xpVyFocalMCWgt,"hsyfp","xpfp_mrad","mcWgt");
  auto h2_xpVypFocalMCWgt= d3.Histo2D(m_xpVypFocalMCWgt,"ypfp_mrad","xpfp_mrad","mcWgt");
  auto h2_yVypFocalMCWgt = d3.Histo2D(m_yVypFocalMCWgt,"ypfp_mrad","hsyfp","mcWgt");
  auto h2_yVxpTarMCWgt   = d3.Histo2D(m_yVxpTarMCWgt,"xptar_mrad","hsytar","mcWgt");
  auto h2_yVypTarMCWgt   = d3.Histo2D(m_yVypTarMCWgt,"yptar_mrad","hsytar","mcWgt");
  auto h2_xpVypTarMCWgt  = d3.Histo2D(m_xpVypTarMCWgt,"yptar_mrad","xptar_mrad","mcWgt");

  auto h2_deltaVypTarMCWgt = d3.Histo2D(m_deltaVypTarMCWgt,"yptar_mrad","delta","mcWgt");

  //XpTar Vs. xfp and xpfp for strange event check
  //auto h2_xpTarVxpFocalMCWgt = d3.Histo2D(m_xpTarVxpFocalMCWgt, "xpfp_mrad", "xptar_mrad", "mcWgt");
  //auto h2_xpTarVxFocalMCWgt = d3.Histo2D(m_xpTarVxFocalMCWgt, "psxfp", "xptar_mrad", "mcWgt");
  //Proj. to SCIN planes
  auto h2_xVyS1XMCWgt = d3.Histo2D(m_xVyS1XMCWgt, "S1X_Y", "S1X_X", "mcWgt");
  auto h2_xVyS1YMCWgt = d3.Histo2D(m_xVyS1YMCWgt, "S1Y_Y", "S1Y_X", "mcWgt");
  auto h2_xVyS2XMCWgt = d3.Histo2D(m_xVyS2XMCWgt, "S2X_Y", "S2X_X", "mcWgt");
  auto h2_xVyS2YMCWgt = d3.Histo2D(m_xVyS2YMCWgt, "S2Y_Y", "S2Y_X", "mcWgt");
  auto h2_dipoleExitMCWgt  = d3.Histo2D(m_dipoleExitMCWgt , "yDipoleExit","xDipoleExit"   ,"mcWgt");
  auto h2_caloEntrMCWgt  = d3.Histo2D(m_caloEntrMCWgt , "yCaloEntr","xCaloEntr"   ,"mcWgt");
  auto h2_ngcerEntrMCWgt  = d3.Histo2D(m_ngcerEntrMCWgt , "yNGCEREntr","xNGCEREntr"   ,"mcWgt");
  auto h2_hgcerEntrMCWgt  = d3.Histo2D(m_hgcerEntrMCWgt , "yHGCEREntr","xHGCEREntr"   ,"mcWgt");


  mcWgtDir = dynamic_cast <TDirectory*> (compFile->Get("mcWgtDir"));
  if(!mcWgtDir) {mcWgtDir = compFile->mkdir("mcWgtDir"); mcWgtDir->cd();}
  else if(mcWgtDir) {compFile->rmdir("mcWgtDir"); mcWgtDir = compFile->mkdir("mcWgtDir"); mcWgtDir->cd();}

  //Write all the histograms to the output file.
  //Define all the histograms
  //h2_xpTarVxpFocalMCWgt->Write();
  //h2_xpTarVxFocalMCWgt->Write();
  h2_xVyS1XMCWgt->Write();
  h2_xVyS1YMCWgt->Write();
  h2_xVyS2XMCWgt->Write();
  h2_xVyS2YMCWgt->Write();
  h2_dipoleExitMCWgt->Write();
  h2_caloEntrMCWgt->Write();
  h2_ngcerEntrMCWgt->Write();
  h2_hgcerEntrMCWgt->Write();
  h2_deltaVypTarMCWgt->Write();
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
  // Weighted 2D MC histos
  h2_xVxpFocalMCWgt->Write();
  h2_xVyFocalMCWgt->Write();
  h2_xVypFocalMCWgt->Write();
  h2_xpVyFocalMCWgt->Write();
  h2_xpVypFocalMCWgt->Write();
  h2_yVypFocalMCWgt->Write();
  h2_yVxpTarMCWgt->Write();
  h2_yVypTarMCWgt->Write();
  h2_xpVypTarMCWgt->Write();

  return 0;
}
