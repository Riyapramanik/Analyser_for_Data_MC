#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TChain.h>

using namespace std;

class PhotonMETSkimmer {
private:
  Int_t nPho;
  vector<float>* pho_pT;
  Float_t met_pT;

  TFile* outputFile;
  TTree* outputTree;

  Long64_t totalEvents;
  Long64_t passedEvents;

public:
  PhotonMETSkimmer() :
    pho_pT(nullptr), outputFile(nullptr), outputTree(nullptr),
    totalEvents(0), passedEvents(0) {}

  ~PhotonMETSkimmer() {
    cleanup();
  }
  
  void cleanup() {
    if (outputFile && outputFile->IsOpen()) {
      outputFile->Close();
      delete outputFile;
      outputFile = nullptr;
    }
  }

  int getValidLeadingPhotonIndex() {
    if (nPho < 1) return -1;
    if (!pho_pT || pho_pT->empty()) return -1;
  
    std::vector<std::pair<float, int>> photonPtIndex;
    for (int i = 0; i < nPho && i < pho_pT->size(); i++) {
        photonPtIndex.push_back(std::make_pair(pho_pT->at(i), i));
    }
    
    std::sort(photonPtIndex.begin(), photonPtIndex.end(), 
              std::greater<std::pair<float, int>>());
    if (photonPtIndex[0].first <= 200.0) return -1;
    
    if (photonPtIndex.size() > 1) {
        if (photonPtIndex[1].first > 200.0) return -1;
    }    
    return photonPtIndex[0].second;
    }

  /* int getValidLeadingPhotonIndex() {
    // Check exactly one photon in event with pT > 200 GeV
    if (nPho != 1) return -1;
    if (!pho_pT || pho_pT->empty()) return -1;
    if (pho_pT->at(0) <= 200.0) return -1;
    
    return 0;
    }*/

  bool passMETCuts() {
    return (met_pT > 200.0);
  }

  void processFileListOptimized(const std::string& fileListPath, const std::string& outputFileName) {
    ifstream fileList(fileListPath);
    TChain* inputChain = new TChain("ggNtuplizer/EventTree");
    
    vector<string> inputFiles;
    string line;
    int validFiles = 0;

    while (getline(fileList, line)) {
      line.erase(0, line.find_first_not_of(" \t\r\n"));
      line.erase(line.find_last_not_of(" \t\r\n") + 1);
      
      if (!line.empty() && line[0] != '#') {
        int result = inputChain->Add(line.c_str());
        validFiles++;
        inputFiles.push_back(line);
      }
    }
    fileList.close();

    pho_pT = nullptr;
    inputChain->SetBranchAddress("nPho", &nPho);
    inputChain->SetBranchAddress("phoPt", &pho_pT);
    inputChain->SetBranchAddress("PuppiMET_pt", &met_pT);

    outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    outputTree = inputChain->CloneTree(0);
    
    Long64_t nentries = inputChain->GetEntries();
    
    for (Long64_t i = 0; i < nentries; i++) {
      inputChain->GetEntry(i);
      totalEvents++;

      int leadingPhotonIndex = getValidLeadingPhotonIndex();
      if (leadingPhotonIndex >= 0 && passMETCuts()) {
        outputTree->Fill();
        passedEvents++;
      }
    }

    TH1F* eventCounts = new TH1F("EventCounts", "Event Counts;Cut;Events", 2, 0, 2);
    eventCounts->SetBinContent(1, totalEvents);
    eventCounts->SetBinContent(2, passedEvents);
    eventCounts->GetXaxis()->SetBinLabel(1, "Total");
    eventCounts->GetXaxis()->SetBinLabel(2, "Passed");

    outputFile->cd();
    outputTree->Write();
    eventCounts->Write();
    outputFile->Close();
    
    delete outputFile;
    outputFile = nullptr;
    delete inputChain;
  }

  void processFileList(const std::string& fileListPath, const std::string& outputFileName) {
    processFileListOptimized(fileListPath, outputFileName);
  }
};

int main(int argc, char* argv[]) {
    string outputFileName = argv[1];
    string fileListPath = argv[2];
    PhotonMETSkimmer skimmer;
    skimmer.processFileList(fileListPath, outputFileName);
    
    skimmer.cleanup();
    exit(0);
        
    return 0;
}
