#pragma once

#define D_CALO_FP 292.64
#define D_EXIT_FP -307.
#define D_S1X_FP 52.1
#define D_S1Y_FP 61.7
#define D_S2X_FP 271.4
#define D_S2Y_FP 282.4
#define D_NGCER_FP -291.700
#define D_HGCER_FP 72.600

using namespace ROOT::RDF;

int weightMC(string requestedKinematic, string target);

//============================================
//        Initialize Physical Constants
//============================================
double sigave = 0.0;

double GetInpVar(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    double var;
    file >> var;
    return var;
}

/*
bool GetKinematicInfoByString(std::fstream& file, const std::string& searchString,
                     std::string& var0, float& var2, float& var3, float& var4) {
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
	      if (!(iss >> var0 >> var2 >> var3 >> var4)) {
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
*/

bool GetTargetInfoByString(std::fstream& file, const std::string& searchString,
			   int& var0, float& var1, float& var2, float& var3, int& var4,
			   float& var5, float& var6, float& var7) {
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
	      if (!(iss >> var0 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7)) {
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

const char* targetInfoFilename = "./target/target_ladder_nps.dat";

//string spectrometer, target;
//float eBeam, eCentral, thetaCentral;


//const Int_t NBINS = 28;
//Double_t edges[NBINS + 1] = {0.35,0.4,0.45,0.50,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95, 1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35, 1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1};

TDirectory *mcWgtDir;

//1D Focal Plane Distributions
  TH1DModel m_xFocalMCWgt ("h_xFocalMCWgt","Weighted Monte-Carlo: X_{fp}; X_{fp} (cm); Number of Entries / 5 mm",100, -40, 40);
  TH1DModel m_xpFocalMCWgt("h_xpFocalMCWgt","Weighted Monte-Carlo: X'_{fp}; X'_{fp}; Number of Entries / 2 mrad",100, -100.0, 100.0);
  TH1DModel m_yFocalMCWgt ("h_yFocalMCWgt","Weighted Monte-Carlo: Y_{fp}; Y_{fp} (cm); Number of Entries / 5 mm",100, -40, 40);
  TH1DModel m_ypFocalMCWgt("h_ypFocalMCWgt","Weighted Monte-Carlo: Y'_{fp}; Y'_{fp}; Number of Entries / 2 mrad",100, -100.0, 100.0);
//1D Target Reconstructed Distributions
  TH1DModel m_yTarMCWgt   ("h_yTarMCWgt","Weighted Monte-Carlo: Y_{tar}; Y_{tar} (cm); Number of Entries / 1 mm",334, -10, 10);
  TH1DModel m_xpTarMCWgt  ("h_xpTarMCWgt","Weighted Monte-Carlo: X'_{tar}; X'_{tar}; Number of Entries / 2 mrad",100, -100.0, 100.0);
  TH1DModel m_ypTarMCWgt  ("h_ypTarMCWgt","Weighted Monte-Carlo: Y'_{tar}; Y'_{tar}; Number of Entries / 2 mrad",100, -100.0, 100.0);
  TH1DModel m_deltaMCWgt  ("h_deltaMCWgt","Weighted Monte-Carlo: #delta; #delta; Number of Entries",120, -30, 30);
//Kinematic Quantities
  TH1DModel m_q2MCWgt     ("h_q2MCWgt","Weighted Monte-Carlo: Q^{2}; Q^{2} (GeV^{2}); Number of Entries / 0.025 GeV^{2}",240, 0.0, 6.0);
  TH1DModel m_w2MCWgt     ("h_w2MCWgt","Weighted Monte-Carlo: W^{2}; W^{2} (GeV^{2}); Number of Entries / 0.050 GeV^{2}",375, -10.0, 20.0);
//TH1DModel m_xbjMCWgt    ("h_xbjMCWgt","Weighted Monte-Carlo: X_{bj}; X_{bj}; Number of Entries / 0.050", NBINS, edges);
TH1DModel m_xbjMCWgt    ("h_xbjMCWgt","Weighted Monte-Carlo: X_{bj}; X_{bj}; Number of Entries / 0.050", 320, 0.0,3.0);
//2D Focal Plane Distributions
  TH2DModel m_xVxpFocalMCWgt ("h2_xVxpFocalMCWgt","Weighted Monte-Carlo: X_{fp} vs. X'_{fp}; X'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",
			       100, -100.0, 100.0, 160, -40, 40);
  TH2DModel m_xVypFocalMCWgt ("h2_xVypFocalMCWgt","Weighted Monte-Carlo: X_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X_{fp} (cm) / 5 mm",
			       60, -60.0, 60.0, 160, -40, 40);
  TH2DModel m_xVyFocalMCWgt  ("h2_xVyFocalMCWgt","Weighted Monte-Carlo: X_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X_{fp} (cm) / 5 mm",
			       160, -40, 40, 160, -40, 40);
  TH2DModel m_xpVyFocalMCWgt ("h2_xpVyFocalMCWgt","Weighted Monte-Carlo: X'_{fp} vs. Y_{fp}; Y_{fp} (cm) / 5 mm; X'_{fp} / 2 mrad",
			       160, -40, 40, 100, -100.0, 100.0);
  TH2DModel m_xpVypFocalMCWgt("h2_xpVypFocalMCWgt","Weighted Monte-Carlo: X'_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; X'_{fp} / 2 mrad",
			       60,-60.0,60.0,100,-100.0, 100.0);
  TH2DModel m_yVypFocalMCWgt ("h2_yVypFocalMCWgt","Weighted Monte-Carlo: Y_{fp} vs. Y'_{fp}; Y'_{fp} / 2 mrad; Y_{fp} (cm) / 5 mm",
			       60,-60.0,60.0,160,-40,40);
//2D Target Distributions
  TH2DModel m_yVxpTarMCWgt   ("h2_yVxpTarMCWgt","Weighted Monte-Carlo: Y_{tar} vs. X'_{tar}; X'_{tar} / 2 mrad; Y_{tar} / 1 mm",
			       200,-100.0,100.0,100,-5,5);
  TH2DModel m_yVypTarMCWgt   ("h2_yVypTarMCWgt","Weighted Monte-Carlo: Y_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; Y_{tar} / 1 mm",
			       200,-100.0,100.0,100,-5,5);
  TH2DModel m_xpVypTarMCWgt  ("h2_xpVypTarMCWgt","Weighted Monte-Carlo: X'_{tar} vs. Y'_{tar}; Y'_{tar} / 2 mrad; X'_{tar} / 2 mrad",
			       60,-60.0,60.0,100,-100.0, 100.0);
//Delta vs. yptar depedence
  TH2DModel m_deltaVypTarMCWgt   ("h2_deltaVypTarMCWgt"   ,"Data: #delta vs. Y'_{Tar}; Y'_{Tar}; #delta",100,-100,100,428, -10,22);
//SCIN locatinos
  TH2DModel m_xVyS1XMCWgt ("h2_xVyS1XMCWgt" ,"MCWgt: X_{S1X} vs. Y_{S1X}; Y_{S1X} (cm) / 5 mm; X_{S1X} (cm) / 5 mm"   ,200,-100.0,100.0,200, -100  , 100  );
  TH2DModel m_xVyS1YMCWgt ("h2_xVyS1YMCWgt" ,"MCWgt: X_{S1Y} vs. Y_{S1Y}; Y_{S1Y} (cm) / 5 mm; X_{S1Y} (cm) / 5 mm"   ,200,-100.0,100.0,200, -100  , 100  );
  TH2DModel m_xVyS2XMCWgt ("h2_xVyS2XMCWgt" ,"MCWgt: X_{S2X} vs. Y_{S2X}; Y_{S2X} (cm) / 5 mm; X_{S2X} (cm) / 5 mm"   ,200,-100.0,100.0,200, -100  , 100  );
  TH2DModel m_xVyS2YMCWgt ("h2_xVyS2YMCWgt" ,"MCWgt: X_{S2Y} vs. Y_{S2Y}; Y_{S2Y} (cm) / 5 mm; X_{S2Y} (cm) / 5 mm"   ,200,-100.0,100.0,200, -100  , 100  );
  TH2DModel m_dipoleExitMCWgt ("h2_dipoleExitMCWgt" ,"Tracks Projected to Dipole Exit; Y_{pos}; X_{pos}"   ,200, -100.0, 100.0,200, -100.0, 100.0);
  TH2DModel m_caloEntrMCWgt ("h2_caloEntrMCWgt" ,"Tracks Projected to Calo Face; Y_{pos}; X_{pos}"          ,200, -100.0, 100.0,200, -100.0, 100.0);
  TH2DModel m_ngcerEntrMCWgt ("h2_ngcerEntrMCWgt" ,"Tracks Projected to NGCER Face; Y_{pos}; X_{pos}" ,200, -100.0, 100.0,200, -100.0, 100.0);
  TH2DModel m_hgcerEntrMCWgt ("h2_hgcerEntrMCWgt" ,"Tracks Projected to HGCER Face; Y_{pos}; X_{pos}" ,200, -100.0, 100.0,200, -100.0, 100.0);
