#include <iostream>
#include <string>
#include <TCutG.h>
#include <TF1.h>
#include <ROOT/RDataFrame.hxx>

string spectrometer, target;
float eBeam, eCentral, thetaCentralRadians, thetaCentral;

string inputFilePattern = "./data/replays/nps_hms_skim_%d_1_-1.root";
string inputFilePattern2 = "./data/replays/nps_hms_skim_%d_1_-1.root_2nd";
//string inputFilePattern = "./data/replays/hms_replay_production_%d_-1.root";
//Pattern is requestedKinematic_runNumber.root
string outputFilePattern = "data/normalized/%s_%d.root";

const char* kinematicFilename = "./kinematics/kinSettingsAll.dat";

bool GetKinematicInfoByString(std::fstream& file, const std::string& searchString,
			      std::string& var1, float& var2, float& var3, float& var4) {
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
	      if (!(iss >> var1 >> var2 >> var3 >> var4)) {
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

bool GetTargetsInKinematic(std::string requestedKinematic, std::set<std::string>& uniqueWords) {
  std::ifstream file(Form("kinematics/%s.dat",requestedKinematic.c_str())); // Replace "data.txt" with your file name
    std::string line;
    //std::set<std::string> uniqueWords;

    if (file.is_open()) {
        while (getline(file, line)) {
            std::istringstream iss(line);
            std::string word;
            int wordCount = 0;

            // Read each word in the line
            while (iss >> word) {
                wordCount++;
                if (wordCount == 2) { // Check if it's the second word
                    uniqueWords.insert(word);
                    break; // Move to the next line after processing the second word
                }
            }
        }
        file.close();

        // Output all unique second words
        std::cout << "Unique targets in this setting:" << std::endl;
        for (const auto& entry : uniqueWords) {
            std::cout << entry << std::endl;
        }
    } else {
        std::cerr << "Unable to open file" << std::endl;
	return false;
    }

    return true;
}

bool GetRunsInKinematic(std::string requestedKinematic, std::vector<int>& runs) {
  std::ifstream file(Form("kinematics/%s.dat",requestedKinematic.c_str())); // Replace "data.txt" with your file name
    std::string line;
    //std::set<std::string> runs;

    if (file.is_open()) {
        while (getline(file, line)) {
            std::istringstream iss(line);
            int word;
            int wordCount = 0;

            // Read each word in the line
            while (iss >> word) {
                wordCount++;
                if (wordCount == 1) { // Check if it's the second word
                    runs.push_back(word);
                    break; // Move to the next line after processing the second word
                }
            }
        }
        file.close();

        // Output all unique second words
        std::cout << "Unique targets in this setting:" << std::endl;
        for (const auto& entry : runs) {
            std::cout << entry << std::endl;
        }
    } else {
        std::cerr << "Unable to open file" << std::endl;
	return false;
    }

    return true;
}


bool GetRunInfoByString(std::fstream& file, std::vector<std::pair<int, std::string>>& data) {
    std::string line;
    // Move back to the beginning of the file
    file.seekg(0, std::ios::beg);
    int num;
    std::string str;    
    std::string loop;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] == '#') {
            continue;
        }
        std::istringstream iss(line);
	if (!(iss >> num >> str  >> loop)) {
	  // Handle error if the line doesn't have enough data
	  std::cerr << "Error: Line format incorrect." << std::endl;
	  return false;
	}
	data.push_back(std::make_pair(num, str));
    }
    return true; // Found the string and extracted the data
}
