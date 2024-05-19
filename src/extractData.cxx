#include <iostream>
#include <string>
#include <TCutG.h>
#include <TF1.h>
#include <ROOT/RDataFrame.hxx>

#include "extractData.h"
#include "processRun.cxx"
//#include "readutils.h"
//#include "subtraction.h"
//#include "directories.h"
//#include "getPionContamination.cxx"
//#include "radCorTbl2.h"

using namespace ROOT::RDF;
using namespace std;

bool extractData(string requestedKinematic) {

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

//Searching for the list of runs for this kinematic (file with first column int (run numbers))
  const char* runlistDirectory = "./kinematics/%s.dat";
  std::vector<std::pair<int, std::string>> runInfo;
  string runlistFilename = Form(runlistDirectory, requestedKinematic.c_str());
  std::cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  std::cout << "                 Run List\n";
  std::cout << "=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:\n";
  std::cout << "Sourcing RunList from:\n" << runlistFilename << endl;
  std::fstream runlistFile(runlistFilename, std::ios::in);
  if (!runlistFile.is_open()) {
    std::cerr << "Error: Could not open file." << std::endl;
    return -1; // Exit with error code if file cannot be opened
  }
  if (!GetRunInfoByString(runlistFile, runInfo)) {
    runlistFile.close();
    return -2; // Exit with a custom error code if string not found or other error
  }

  int runNumber; std::string target;
  //Process each run in the runlist and save the output to a unique TFile
  for (const auto& run : runInfo) { 
    runNumber = run.first;
    target = run.second;
    //cout << runNumber << " " << target << endl;
    //continue;
    //Check the root file is valid for a given run number.

//If the run number is not valid, it could be a junk run, wrong path, wrong filename.
//More likely the file is not in /cache/.  Stop analysis and get with jcache
//Make sure T Tree exists in the input file.
//This is backwards from how the input files should be searched.  Look for _2nd first.  
    string inputFilename = Form(inputFilePattern2.c_str(), runNumber);
    TFile* inputFile = TFile::Open(inputFilename.c_str(), "READ");  //Read, do not modify
    if (!inputFile || inputFile->IsZombie()) {
      std::cerr << "Error opening file.\nChecking for original .root file..." << inputFilename << std::endl;
      //return -1;
    }
    inputFilename = Form(inputFilePattern.c_str(), runNumber);
    inputFile = TFile::Open(inputFilename.c_str(), "READ");  //Read, do not modify
    if (!inputFile || inputFile->IsZombie()) {
      std::cerr << "Error opening file " << inputFilename << std::endl;
      return -1;
    }
    if (!inputFile->Get("T")) {
      std::cerr << "Error, no TTree named T found in " << inputFilename << std::endl;
      return -5;
    }

//Allow the extractData method to check each outputFilename is valid.
    string outputFilename = Form(outputFilePattern.c_str(),requestedKinematic.c_str(), runNumber);
    TFile* outputFile;
    
    outputFile = new TFile(outputFilename.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
      std::cerr << "Error opening file " << outputFilename << std::endl;
      return -1;
    }
    
    if (!processRun(runNumber, inputFile, outputFile, requestedKinematic, 
	  target, spectrometer, eBeam, eCentral, thetaCentralRadians)) {
      std::cerr << "Error processing run " << runNumber << std::endl;
      return -4; // Or handle error as needed
    }
    
//Method to extract scalar normalizations and errors, save in a TTree
    if(! LoadDataIntoROOT("data_norm.csv",outputFilename.c_str(), runNumber)) {
      //return 0;
    }

    NormalizeHistograms(outputFilename.c_str(), runNumber);

  }
  
  return 1;
}
