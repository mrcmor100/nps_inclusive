#include "processRun.h"

using namespace ROOT::RDF;
using namespace std;

//Add a README / HEADER here to explain what is happening.

bool processRun(int run, TFile* inputFile, TFile* outputFile, string requestedKinematic, 
		string target, string spectrometer, float eBeam, float eCentral, float thetaCentralRadians) {

  auto calcCalEff = [=](double deltaData) {
    //TODO: No calorimeter efficiency has been measured for HMS (to my knowledge).
    //SHMS: return 0.9984 - TMath::Exp(-1.98 - 0.2356 * TMath::Power(eCentral * (1 + deltaData / 100), 3));
    return 1.0;
  };
  auto rad2mrad   = [](double pFocalData) {
    return pFocalData*1000.0;
  };
  auto calc_q2 = [=] (double th, double ph, double dp) {
    double hse = eCentral*(1. + dp/100.);
    double nu = eBeam - hse;
    double hstheta = acos(cos(thetaCentralRadians + ph)*cos(th));
    double sn2 = sin(hstheta/2.)*sin(hstheta/2.);
    return 4.*hse*eBeam*sn2;
  };
  auto calc_w2 = [=] (double q2_me, double dp) {
    double hse = eCentral*(1. + dp/100.);
    double nu = eBeam - hse;
    return mp2 + 2.*mp*nu-q2_me;
  };
  auto calc_xbj = [=] (double q2_me, double dp) {
    double hse = eCentral*(1. + dp/100.);
    double nu = eBeam - hse;
    return q2_me / (2*mp*nu);
  };

  //Open TFile here.
  //Check it exists, if not error.
  ROOT::EnableImplicitMT(8);
  //std::cout << ROOT::GetImplicitMTPoolSize() << endl
  TTree* T = (TTree*) inputFile->Get("T");
  ROOT::RDataFrame d(*T);
  //Filter all the crap events BEFORE making histograms and new variables
  
  auto nEntries = d.Count();
  std::cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  std::cout << "                  Total Number of Entries \n";
  std::cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  std::cout << "The total number of Entries is: " << *nEntries << endl;
  std::cout << endl << endl;

  auto d2 = d.Filter("H.cer.npeSum > 0.7")
    .Filter("H.cal.etottracknorm > 0.7")
    .Filter("H.gtr.dp > -8. && H.gtr.dp < 8.")
    .Filter("abs(H.gtr.ph) < 0.050")
    .Filter("abs(H.gtr.th) < 0.22");
    //.Filter("abs(H.gtr.y) < 4.")
    //.Filter("P.bcm.bcm4c.AvgCurrent > 5")
    //.Filter("P.dc.InsideDipoleExit==1")

  //Define new variables of interest.
  d2 = d2.Define("P_gtr_xp_mrad" ,rad2mrad, {"H.gtr.th"}  )
    .Define("P_gtr_yp_mrad" ,rad2mrad, {"H.gtr.ph"}  )
    .Define("P_dc_xpfp_mrad",rad2mrad, {"H.dc.xp_fp"})
    .Define("P_dc_ypfp_mrad",rad2mrad, {"H.dc.yp_fp"})
    .Define("q2_me", calc_q2, {"H.gtr.th","H.gtr.ph","H.gtr.dp"})
    .Define("w2_me", calc_w2, {"q2_me","H.gtr.dp"})
    .Define("w_me", "sqrt(abs(w2_me))")
    .Define("xbj_me", calc_xbj, {"q2_me","H.gtr.dp"});

  auto nEvents = d2.Count();
  std::cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  std::cout << "             Good Electron Events\n";
  std::cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  std::cout << "Number of good electron events: " << *nEvents << 
    endl << endl << endl;

  auto h_calEttnData  = d2.Histo1D(m_calEttnData  ,"H.cal.etottracknorm");
  auto h_xFocalData  = d2.Histo1D(m_xFocalData  ,"H.dc.x_fp");
  auto h_xpFocalData = d2.Histo1D(m_xpFocalData ,"P_dc_xpfp_mrad");
  auto h_yFocalData  = d2.Histo1D(m_yFocalData  ,"H.dc.y_fp");
  auto h_ypFocalData = d2.Histo1D(m_ypFocalData ,"P_dc_ypfp_mrad");
  auto h_yTarData    = d2.Histo1D(m_yTarData , "H.gtr.y");
  auto h_xpTarData   = d2.Histo1D(m_xpTarData, "P_gtr_xp_mrad");
  auto h_ypTarData   = d2.Histo1D(m_ypTarData, "P_gtr_yp_mrad");
  auto h_deltaData   = d2.Histo1D(m_deltaData, "H.gtr.dp");
  auto h_xbjData     = d2.Histo1D(m_xbjData  , "H.kin.x_bj");
  auto h_xbjDataMe   = d2.Histo1D(m_xbjDataMe, "xbj_me");
  auto h_w2DataMe    = d2.Histo1D(m_w2DataMe,  "w2_me");
  auto h_q2DataMe    = d2.Histo1D(m_q2DataMe,  "q2_me");
  auto h_q2Data      = d2.Histo1D(m_q2Data   , "H.kin.Q2");
  auto h_w2Data      = d2.Histo1D(m_w2Data   , "H.kin.W2");
  auto h_wDataMe     = d2.Histo1D(m_wDataMe  , "w_me");
  auto h_wData       = d2.Histo1D(m_wData    , "H.kin.W");

  // Data 2D histos
  auto h2_xVxpFocalData  = d2.Histo2D(m_xVxpFocalData , "P_dc_xpfp_mrad","H.dc.x_fp");
  auto h2_xVyFocalData   = d2.Histo2D(m_xVyFocalData  , "H.dc.y_fp","H.dc.x_fp");
  auto h2_xVypFocalData  = d2.Histo2D(m_xVypFocalData , "P_dc_ypfp_mrad","H.dc.x_fp");
  auto h2_xpVyFocalData  = d2.Histo2D(m_xpVyFocalData , "H.dc.y_fp","P_dc_xpfp_mrad");
  auto h2_xpVypFocalData = d2.Histo2D(m_xpVypFocalData, "P_dc_ypfp_mrad","P_dc_xpfp_mrad");
  auto h2_yVypFocalData  = d2.Histo2D(m_yVypFocalData , "P_dc_ypfp_mrad","H.dc.y_fp");
  auto h2_yVxpTarData    = d2.Histo2D(m_yVxpTarData   , "P_gtr_xp_mrad","H.gtr.y");
  auto h2_yVypTarData    = d2.Histo2D(m_yVypTarData   , "P_gtr_yp_mrad","H.gtr.y");
  auto h2_xpVypTarData   = d2.Histo2D(m_xpVypTarData  , "P_gtr_yp_mrad","P_gtr_xp_mrad");

  // Create data directory and descend into it
  dataDir = dynamic_cast <TDirectory*> (outputFile->Get("dataHistograms"));
  if(!dataDir) {
    dataDir = outputFile->mkdir("dataHistograms");
    dataDir->cd();
  } else if(dataDir) {
    outputFile->rmdir("dataHistograms");
    dataDir = outputFile->mkdir("dataHistograms");
    dataDir->cd();
  }

  h_calEttnData->Write(); 
  h_xFocalData->Write();
  h_xpFocalData->Write();
  h_yFocalData->Write();
  h_ypFocalData->Write();
  h_yTarData->Write();
  h_xpTarData->Write();
  h_ypTarData->Write();
  h_deltaData->Write();
  //Kinematic Variables
  h_xbjData->Write();
  h_xbjDataMe->Write();
  h_q2Data->Write();
  h_q2DataMe->Write();
  h_wData->Write();
  h_wDataMe->Write();
  h_w2Data->Write();
  h_w2DataMe->Write();
  //2D focal plane variables
  h2_xVxpFocalData->Write();
  h2_xVyFocalData->Write();
  h2_xVypFocalData->Write();
  h2_xpVyFocalData->Write();
  h2_xpVypFocalData->Write();
  h2_yVypFocalData->Write();
  h2_yVxpTarData->Write();
  h2_yVypTarData->Write();
  h2_xpVypTarData->Write();

  outputFile->Close();
  return true;
}

bool addWeightedHistograms(int run, TFile* inputFile, TFile* outputFile, string requestedKinematic, 
			   string target, string spectrometer, float eBeam, float eCentral, float thetaCentralRadians) {

  // Access the tree
  //TTree* tree = (TTree*)file->Get("scalarTree");
  //if (!tree) {
  //std::cerr << "Tree not found!" << std::endl;
  //file->Close();
  //return;
  //}


  return true;
}

/*

// Function to read scalar values from an existing ROOT file into the variables.
bool readScalars(const char* filename) {
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    TTree* tree;
    file->GetObject("ScalarsTree", tree);
    if (!tree) {
        std::cerr << "Tree not found in file!" << std::endl;
        file->Close();
        return false;
    }

    // Set branch addresses
    tree->SetBranchAddress("ScalarValue1", &scalar1);
    tree->SetBranchAddress("ScalarValue2", &scalar2);
    tree->SetBranchAddress("ScalarValue3", &scalar3);

    // Read the first entry (assuming only one entry is present)
    tree->GetEntry(0);

    file->Close();
    return true;
}
*/
/*
bool LoadDataIntoROOT(const char* csvFilePath, const char* rootFilePath, int runNumber) {
    // Open the CSV file
    std::ifstream csvFile(csvFilePath);
    std::string line, cell;
    
    // Open the ROOT file
    TFile *file = new TFile(rootFilePath, "UPDATE");
    TTree *tree = new TTree("ScalarTree", "Tree with scalar data from CSV");

    // Check if the file is open
    if (!csvFile.is_open()) {
        std::cerr << "Failed to open CSV file!" << std::endl;
        return false;
    }
    
    // Read the header from the CSV
    std::getline(csvFile, line);
    std::istringstream headerStream(line);
    std::vector<std::string> headers;
    while (std::getline(headerStream, cell, ',')) {
        headers.push_back(cell);
    }
    
    // Prepare to hold the data
    std::vector<float> data(headers.size(), 0.0);
    std::vector<TBranch*> branches;

    // Create branches for each column except the first (run number)
    for (size_t i = 1; i < headers.size(); ++i) {
        branches.push_back(tree->Branch(headers[i].c_str(), &data[i]));
    }
    
    bool foundRun = false;
    // Read the data from the CSV
    while (std::getline(csvFile, line)) {
        std::istringstream lineStream(line);
        int currentRunNumber;
        lineStream >> currentRunNumber;
        // Check if this line corresponds to the specified run number
        if (currentRunNumber == runNumber) {
	  foundRun=true;
            size_t index = 0;
            while (std::getline(lineStream, cell, ',')) {
                data[index] = std::stof(cell);
                index++;
            }
            
            // Fill the tree
            tree->Fill();
            break; // Assuming only one record per run number
        }
	if(!foundRun)
	  return false;
    }

    // Write the tree to file
    tree->Write();
    file->Close();

    // Clean up
    delete file;
    csvFile.close();

    return true;
}
*/
bool LoadDataIntoROOT(const char* csvFilePath, const char* rootFilePath, int runNumber) {
    // Open the CSV file
    std::ifstream csvFile(csvFilePath);
    std::string line, cell;

    // Open the ROOT file
    TFile *file = new TFile(rootFilePath, "UPDATE");
    TTree *tree = new TTree("ScalarTree", "Tree with scalar data from CSV");

    // Check if the file is open
    if (!csvFile.is_open()) {
        std::cerr << "Failed to open CSV file!" << std::endl;
        return false;
    }

    // Read the header from the CSV
    std::getline(csvFile, line);
    std::istringstream headerStream(line);
    std::vector<std::string> headers;
    while (std::getline(headerStream, cell, ',')) {
        headers.push_back(cell);
    }

    // Prepare to hold the data
    std::vector<float> data(headers.size(), 0.0);
    std::vector<TBranch*> branches;

    // Create branches for each column except the first (run number)
    for (size_t i = 1; i < headers.size(); ++i) {
        branches.push_back(tree->Branch(headers[i].c_str(), &data[i]));
    }

    bool foundRun = false;

    // Read the data from the CSV
    while (std::getline(csvFile, line)) {
        std::istringstream lineStream(line);
        std::getline(lineStream, cell, ',');
        int currentRunNumber = std::stoi(cell);

        // Check if this line corresponds to the specified run number
        if (currentRunNumber == runNumber) {
            foundRun = true;
            size_t index = 1;  // Start from 1 because index 0 is run number
            while (std::getline(lineStream, cell, ',')) {
                data[index] = std::stof(cell);
                index++;
            }

            // Fill the tree
            tree->Fill();
            break; // Assuming only one record per run number
        }
    }

    if (!foundRun) {
        std::cerr << "Run number " << runNumber << " not found in CSV file!" << std::endl;
        return false;
    }

    // Write the tree to file
    tree->Write();
    file->Close();

    // Clean up
    delete file;
    csvFile.close();

    return true;
}

void NormalizeHistograms(const char* rootFilePath, int runNumber) {
    // Open the ROOT file
    TFile *file = new TFile(rootFilePath, "UPDATE");
    if (!file->IsOpen()) {
        std::cerr << "Failed to open ROOT file!" << std::endl;
        return;
    }

    // Load the tree and set up branches
    TTree *tree = (TTree*)file->Get("ScalarTree");
    if (!tree) {
        std::cerr << "Failed to find ScalarTree in ROOT file!" << std::endl;
        return;
    }

    float BCM4A, Effic, Ps3Time, Ps4Time;
    float Ps3Factor, Ps4Factor;
    tree->SetBranchAddress("BCM4AQ", &BCM4A);
    tree->SetBranchAddress("trkEff", &Effic);
    tree->SetBranchAddress("Ps3HMSCLT", &Ps3Time);
    tree->SetBranchAddress("Ps4HMSCLT", &Ps4Time);
    tree->SetBranchAddress("Ps3", &Ps3Factor);
    tree->SetBranchAddress("Ps4", &Ps4Factor);

    tree->GetEntry(0);

    // Check the Ps factors
    float liveTime = 0.0;
    if ((Ps3Factor > 0 && Ps4Factor > 0) || (Ps3Factor <= 0 && Ps4Factor <= 0)) {
        std::cerr << "Error with prescale factors!" << std::endl;
        return;
    }

    int PsFactor = -10;
    if (Ps3Factor > 0) {
        liveTime = Ps3Time;
	PsFactor = (int)Ps3Factor;
    } else if (Ps4Factor > 0) {
        liveTime = Ps4Time;
	PsFactor = (int)Ps4Factor;
    }

    // Calculate the normalization factor
    float normFactor = 1.0 / (BCM4A * Effic * liveTime / PsFactor);

    // Normalize histograms
    TDirectory *histDir = (TDirectory*)file->Get("dataHistograms");
    TDirectory *normDir;// = file->mkdir("dataHistogramsNorm");
    normDir = dynamic_cast <TDirectory*> (file->Get("dataHistogramsNorm"));
    if(!normDir) {
      normDir = file->mkdir("dataHistogramsNorm");
      normDir->cd();
    } else if(dataDir) {
      file->rmdir("dataHistogramsNorm");
      normDir = file->mkdir("dataHistogramsNorm");
      normDir->cd();
    }
    normDir->cd();

    TIter next(histDir->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();
        if (TH1 *h1 = dynamic_cast<TH1*>(obj)) {
            TH1 *h1Clone = (TH1*)h1->Clone();
            h1Clone->Scale(normFactor);
            h1Clone->Write();
        } else if (TH2 *h2 = dynamic_cast<TH2*>(obj)) {
            TH2 *h2Clone = (TH2*)h2->Clone();
            h2Clone->Scale(normFactor);
            h2Clone->Write();
        }
    }

    // Save changes and close file
    //file->Write();
    file->Close();
}
