#include "extractData.h"
#include "dataToMonteCarlo.cxx"

vector<int> dummy_runs;

bool comparisons(std::string requestedKinematic) {

  const char* runlistDirectory = "./kinematics/%s.dat"; //DUPLICATE

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

  string loop;
  int runNumber; std::string target;
  //Process each run in the runlist and save the output to a unique TFile
  for (const auto& run : runInfo) { 
    runNumber = run.first;
    target = run.second;
      
    //Accumulate dummy runs for combining afterwards main loop
    //For now, just use dummy run closest in run number
    cout << target << endl;
    if(target == "dummy") {
      cout << "Adding dummy run to dummy list\n";
      dummy_runs.push_back(runNumber);
    }
  }

  if (!GetRunInfoByString(runlistFile, runInfo)) {
    runlistFile.close();
    return -2; // Exit with a custom error code if string not found or other error
  }
  //Process each run in the runlist and save the output to a unique TFile

  for (const auto& run : runInfo) { 
    runNumber = run.first;
    target = run.second;
    if(target == "ld2") {
      loop = "loop1";
    } else if (target == "lh2") {
      loop = "loop2";
    } else {
      loop = "na";
    }

    int closest = 0;
    if (!dummy_runs.empty()) {
      closest = *std::min_element(dummy_runs.begin(), dummy_runs.end(), [&runNumber](int a, int b) {
	  return std::abs(a - runNumber) < std::abs(b - runNumber);
	});
      cout << "Closest Dummy run to: " << runNumber << " is: " << closest << endl;
    }
    //dummy_runs.clear();

    cout << requestedKinematic << " " << runNumber << " " << target << " " <<loop << " " <<closest << endl;
    dataToMonteCarlo(requestedKinematic, runNumber, target, loop, closest);
  }
  return true;
}
