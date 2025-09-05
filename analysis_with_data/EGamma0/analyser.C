#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TMath.h>
#include <TROOT.h>
#include "Utilities.h"
#include "algorithm"
using namespace std;

class PhotonAnalyser {

private:
  //Photon variable
  Int_t nPho;
  vector<float>* phoPt;
  vector<float>* phoPhi;
  vector<float>* phoSCEta;
  vector<float>* phoMIPTotEnergy;
  vector<float>* phoSigmaIEtaIEtaFull5x5;
  vector<float>* phoSigmaIPhiIPhiFull5x5;
  vector<float>* phoSeedTime;
  vector<int>* phohasPixelSeed;
  vector<float>* phoSCEtaWidth;
  vector<float>* phoHaloTaggerMVAVal;
  vector<UShort_t>* phoIDbit;

  //Jet variables
  UShort_t nAK4PUPPIJet;
  vector<float>* AK4PUPPIJet_Pt;
  vector<float>* AK4PUPPIJet_Eta;
  vector<float>* AK4PUPPIJet_Phi;
  vector<float>* AK4PUPPIJet_jetID;

  // MET variables
  Float_t PuppiMET_pt;
  Float_t PuppiMET_phi;

  // HLT trigger
  UShort_t HLTriggerWord0;

  // MET filter variables
  Bool_t Flag_goodVertices;
  Bool_t Flag_globalSuperTightHalo2016Filter;
  Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter;
  Bool_t Flag_BadPFMuonFilter;
  Bool_t Flag_BadPFMuonDzFilter;
  Bool_t Flag_hfNoisyHitsFilter;
  Bool_t Flag_eeBadScFilter;

  TFile* outputFile;

  // Histograms
  TH1F* h_pT_pho_EB;
  TH1F* h_pT_pho_EE;
  TH1F* h_eta_pho_EB;
  TH1F* h_eta_pho_EE;
  TH1F* h_phi_pho_EB;
  TH1F* h_phi_pho_EE;
  TH1F* h_pt_Met_EB;
  TH1F* h_phi_Met_EB;
  TH1F* h_pt_Met_EE;
  TH1F* h_phi_Met_EE;
  TH1F* h_jet_deltaPhi_MET;
  TH1F* h_deltaR_pho_jet;
  TH1F* h_cutflow;
  TH1F* h_pho_etawidth_EB;
  TH1F* h_pho_etawidth_EE;
    
  // Counters
  int total_events;
  int hlt_pass;
  int met_filters_pass;
  int pt_pass;
  int pixel_seed_veto_pass;
  int beamhalo_pass;
  int tight_id_pass;
  int met_cuts_pass;
  int jet_veto_pass;
  int barrel_pass;
  int endcap_pass;
  int barrel_MIP_pass;
  int barrel_passSigmaIetaIeta;
  int barrel_passSigmaIphiIphi;
  int barrel_passEtaWidth;
  int barrel_passTiming;
  int barrel_all_cuts;

public:
  PhotonAnalyser() :
    phoPt(nullptr), phoPhi(nullptr), phoSCEta(nullptr), phoMIPTotEnergy(nullptr),
    phoSigmaIEtaIEtaFull5x5(nullptr), phoSigmaIPhiIPhiFull5x5(nullptr), phoSeedTime(nullptr),
    phohasPixelSeed(nullptr), phoSCEtaWidth(nullptr), phoIDbit(nullptr),AK4PUPPIJet_Pt(nullptr),
    AK4PUPPIJet_Eta(nullptr), AK4PUPPIJet_Phi(nullptr), AK4PUPPIJet_jetID(nullptr),
    outputFile(nullptr),
    total_events(0), hlt_pass(0), met_filters_pass(0), pt_pass(0), 
    pixel_seed_veto_pass(0), met_cuts_pass(0), jet_veto_pass(0),
    barrel_pass(0), endcap_pass(0), barrel_all_cuts(0),beamhalo_pass(0), tight_id_pass(0), barrel_MIP_pass(0),barrel_passSigmaIetaIeta(0),barrel_passSigmaIphiIphi(0), barrel_passEtaWidth(0),barrel_passTiming(0) {}
  
  ~PhotonAnalyser() {
    if (outputFile && outputFile->IsOpen()) {
      outputFile->Close();
      delete outputFile;
    }
  }

  void setupInputBranches(TTree* tree) {
    
    // Reset pointers
    phoPt = nullptr;
    phoPhi = nullptr;
    phoSCEta = nullptr;
    phoMIPTotEnergy = nullptr;
    phoSigmaIEtaIEtaFull5x5 = nullptr;
    phoSigmaIPhiIPhiFull5x5 = nullptr;
    phoSeedTime = nullptr;
    phohasPixelSeed = nullptr;
    phoSCEtaWidth = nullptr;
    phoHaloTaggerMVAVal = nullptr;
    phoIDbit = nullptr;
    AK4PUPPIJet_Pt = nullptr;
    AK4PUPPIJet_Eta = nullptr;
    AK4PUPPIJet_Phi = nullptr;
    AK4PUPPIJet_jetID = nullptr;
    
    // Set branch addresses
    tree->SetBranchAddress("nPho", &nPho);
    tree->SetBranchAddress("phoPt", &phoPt);
    tree->SetBranchAddress("phoPhi", &phoPhi);
    tree->SetBranchAddress("phoSCEta", &phoSCEta);
    tree->SetBranchAddress("phoMIPTotEnergy", &phoMIPTotEnergy);
    tree->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5);
    tree->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5);
    tree->SetBranchAddress("phoSeedTime", &phoSeedTime);
    tree->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed);
    tree->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth);
    tree->SetBranchAddress("phoHaloTaggerMVAVal", &phoHaloTaggerMVAVal);
    tree->SetBranchAddress("phoIDbit", &phoIDbit);
    
    // Jet branches
    tree->SetBranchAddress("nAK4PUPPIJet", &nAK4PUPPIJet);
    tree->SetBranchAddress("AK4PUPPIJet_Pt", &AK4PUPPIJet_Pt);
    tree->SetBranchAddress("AK4PUPPIJet_Eta", &AK4PUPPIJet_Eta);
    tree->SetBranchAddress("AK4PUPPIJet_Phi", &AK4PUPPIJet_Phi);
    tree->SetBranchAddress("AK4PUPPIJet_jetID", &AK4PUPPIJet_jetID);
          
    // MET branches
    tree->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt);
    tree->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi);
    
    // HLT trigger
    tree->SetBranchAddress("HLTriggerWord0", &HLTriggerWord0);
          
    // MET filter branches
    tree->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices);
    tree->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter);
    tree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
    tree->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter);
    tree->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter);
    tree->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter);
    tree->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter);
  }

// Function to get the leading photon index
  int getLeadingPhotonIndex() {
    if (nPho < 1) return -1;
    if (!phoPt || phoPt->empty()) return -1;
    std::vector<std::pair<float, int>> photonPtIndex;
    for (int i = 0; i < nPho && i < phoPt->size(); i++) {
        photonPtIndex.push_back(std::make_pair(phoPt->at(i), i));
    }
    std::sort(photonPtIndex.begin(), photonPtIndex.end(), 
              std::greater<std::pair<float, int>>());
    return photonPtIndex[0].second;  
  }

  void setupOutputHistograms(const std::string& outputFileName) {
    outputFile = new TFile(outputFileName.c_str(), "RECREATE");

    const int nBins_pT = 6;
    Double_t bins_pT[7] = {225, 275, 350, 450, 600, 800, 1500};
    const int nBins_met = 5;
    Double_t bins_met[6] = {225, 275, 350, 450, 600, 1000};
    
    // Create histograms
    h_pT_pho_EB = new TH1F("h_pT_EB", "Photon p_{T} - Barrel;p_{T} [GeV];Events", nBins_pT, bins_pT);
    h_pT_pho_EE = new TH1F("h_pT_EE", "Photon p_{T} - Endcap;p_{T} [GeV];Events", nBins_pT, bins_pT);
    h_eta_pho_EB = new TH1F("h_eta_EB", "Photon #eta - Barrel;#eta;Events", 15, -2.0, 2.0);
    h_eta_pho_EE = new TH1F("h_eta_EE", "Photon #eta - Endcap;#eta;Events", 15, -3.5, 3.5);
    h_phi_pho_EB = new TH1F("h_phi_EB", "Photon #phi - Barrel;#phi;Events", 15, -3.5, 3.5);
    h_phi_pho_EE = new TH1F("h_phi_EE", "Photon #phi - Endcap;#phi;Events", 15, -3.5, 3.5);
    
    h_pt_Met_EB = new TH1F("h_MET_pt_EB", "MET p_{T} - Barrel;p_{T} [GeV];Events", nBins_met,bins_met);
    h_phi_Met_EB = new TH1F("h_MET_phi_EB", "MET #phi - Barrel;#phi;Events", 15, -3.5, 3.5);
    h_pt_Met_EE = new TH1F("h_MET_pt_EE", "MET p_{T} - Endcap;p_{T} [GeV];Events", nBins_met,bins_met);
    h_phi_Met_EE = new TH1F("h_MET_phi_EE", "MET #phi - Endcap;#phi;Events", 15, -3.5, 3.5);
    
    h_jet_deltaPhi_MET = new TH1F("h_jet_deltaPhi_MET", "#Delta#phi(jet, MET);#Delta#phi;Events", 10, 0, 3.2);
    h_deltaR_pho_jet = new TH1F("h_deltaR_pho_jet", "#DeltaR(#gamma, jet);#DeltaR;Events", 10, 0, 5);

    h_pho_etawidth_EB = new TH1F("h_pho_etawidth_EB","photon etawidth EB",20,0,0.1);
    h_pho_etawidth_EE =	new TH1F("h_pho_etawidth_EE","photon etawidth EE",20,0,0.1);
    
    h_cutflow = new TH1F("h_cutflow", "Cut Flow;Cut;Events", 17, 0, 17);
  }

  void processEvents(TTree* tree) {
    Long64_t nentries = tree->GetEntries();   
    for (Long64_t i = 0; i < nentries; i++) {
      tree->GetEntry(i);
      total_events++;
      
      // HLT trigger
      bool HLT_Photon200 = (HLTriggerWord0 & (1ULL << 11)) != 0;
      if (!HLT_Photon200) continue;
      hlt_pass++;
      
      // MET filters
      if (!Flag_goodVertices) continue;
      if (!Flag_globalSuperTightHalo2016Filter) continue;
      if (!Flag_EcalDeadCellTriggerPrimitiveFilter) continue;
      if (!Flag_BadPFMuonFilter) continue;
      if (!Flag_BadPFMuonDzFilter) continue;
      if (!Flag_hfNoisyHitsFilter) continue;
      if (!Flag_eeBadScFilter) continue;
      met_filters_pass++;

      // Photon loop
      int leadingPhotonIndex = getLeadingPhotonIndex();
      int iPho = leadingPhotonIndex;
        if (phoPt->at(iPho) <= 225.0) continue;
        pt_pass++;
	
        if (phohasPixelSeed->at(iPho) != 0) continue;
        pixel_seed_veto_pass++;

	if (phoHaloTaggerMVAVal->at(iPho) <= 0.995) continue;
	beamhalo_pass++;

	if(!(((phoIDbit->at(iPho) >> 2) & 1) == 1)) continue;
        tight_id_pass++;

        float pT_over_MET = phoPt->at(iPho) / PuppiMET_pt;
        if (pT_over_MET >= 1.4) continue;

	Double_t deltaPhi_pho_MET = fabs(deltaPhi(phoPhi->at(iPho), PuppiMET_phi));
        if (deltaPhi_pho_MET <= 2.0) continue;
        met_cuts_pass++;
              
        // JET VETO
        bool JetVeto = false;
                  
        for (int iJet = 0; iJet < nAK4PUPPIJet; iJet++) {
	  if (AK4PUPPIJet_Pt->at(iJet) <= 30.0) continue;
	  if (fabs(AK4PUPPIJet_Eta->at(iJet)) >= 5.0) continue;
	  if (AK4PUPPIJet_jetID->at(iJet) < 0.5) continue;
	  
	  Double_t deltaR_pho_jet = deltaR(phoSCEta->at(iPho), phoPhi->at(iPho), 
	                                   AK4PUPPIJet_Eta->at(iJet), AK4PUPPIJet_Phi->at(iJet));
	  
	  h_deltaR_pho_jet->Fill(deltaR_pho_jet);
	  
	  if (deltaR_pho_jet <= 0.5) continue;
	  Double_t deltaPhi_jet_MET = fabs(deltaPhi(AK4PUPPIJet_Phi->at(iJet), PuppiMET_phi));
	  h_jet_deltaPhi_MET->Fill(deltaPhi_jet_MET);
	  
	  if (deltaPhi_jet_MET < 0.5) {
	    JetVeto = true;
	    break;
	  }
        }
	
        if (JetVeto) continue;
        jet_veto_pass++;	
	
        float abs_eta_SC = fabs(phoSCEta->at(iPho));
        bool isBarrel = (abs_eta_SC < 1.4442);
        bool isEndcap = (abs_eta_SC > 1.566 && abs_eta_SC < 2.5);
	
        if (isBarrel) {
	  barrel_pass++;
	  bool passMIP = (phoMIPTotEnergy->at(iPho) < 4.9);
	  bool passSigmaIetaIeta = (phoSigmaIEtaIEtaFull5x5->at(iPho) > 0.001);
	  bool passSigmaIphiIphi = (phoSigmaIPhiIPhiFull5x5->at(iPho) > 0.001);
	  bool passEtaWidth = (phoSCEtaWidth->at(iPho) > 0.01);
	  bool passTiming = (fabs(phoSeedTime->at(iPho)) < 3.0);
	  if (passMIP){barrel_MIP_pass++;}
	  if(passSigmaIetaIeta){barrel_passSigmaIetaIeta++;}
	  if(passSigmaIphiIphi){barrel_passSigmaIphiIphi++;}
	  if(passEtaWidth){barrel_passEtaWidth++;}
	  if(passTiming){barrel_passTiming++;}
	    
	  // if (passMIP && passSigmaIetaIeta && passSigmaIphiIphi && passEtaWidth && passTiming) {
	  if (passMIP && passSigmaIetaIeta && passSigmaIphiIphi && passTiming) {
	    barrel_all_cuts++;
	    h_pT_pho_EB->Fill(phoPt->at(iPho));
	    h_eta_pho_EB->Fill(phoSCEta->at(iPho));
	    h_phi_pho_EB->Fill(phoPhi->at(iPho));
	    
	    h_pt_Met_EB->Fill(PuppiMET_pt);
	    h_phi_Met_EB->Fill(PuppiMET_phi);

	    h_pho_etawidth_EB->Fill(phoSCEtaWidth->at(iPho));
	    
	  }
        } else if (isEndcap) {
	  endcap_pass++;
	  h_pT_pho_EE->Fill(phoPt->at(iPho));
	  h_eta_pho_EE->Fill(phoSCEta->at(iPho));
	  h_phi_pho_EE->Fill(phoPhi->at(iPho));
	  h_pt_Met_EE->Fill(PuppiMET_pt);
	  h_phi_Met_EE->Fill(PuppiMET_phi);
	  h_pho_etawidth_EE->Fill(phoSCEtaWidth->at(iPho));
        }
    }
  }
  
  void saveHistograms() {
    outputFile->cd();
    
    // Fill cut flow histogram
    h_cutflow->SetBinContent(1, total_events);
    h_cutflow->SetBinContent(2, hlt_pass);
    h_cutflow->SetBinContent(3, met_filters_pass);
    h_cutflow->SetBinContent(4, pt_pass);
    h_cutflow->SetBinContent(5, pixel_seed_veto_pass);
    h_cutflow->SetBinContent(6, beamhalo_pass);
    h_cutflow->SetBinContent(7, tight_id_pass);
    h_cutflow->SetBinContent(8, met_cuts_pass);
    h_cutflow->SetBinContent(9, jet_veto_pass);
    h_cutflow->SetBinContent(10, barrel_pass);
    h_cutflow->SetBinContent(11, barrel_MIP_pass);
    h_cutflow->SetBinContent(12, barrel_passSigmaIetaIeta);
    h_cutflow->SetBinContent(13, barrel_passSigmaIphiIphi);
    h_cutflow->SetBinContent(14, barrel_passEtaWidth);
    h_cutflow->SetBinContent(15, barrel_passTiming);
    h_cutflow->SetBinContent(16, barrel_all_cuts);
    h_cutflow->SetBinContent(17, endcap_pass);
    

    h_cutflow->GetXaxis()->SetBinLabel(1, "total_events");
    h_cutflow->GetXaxis()->SetBinLabel(2, "hlt_pass");
    h_cutflow->GetXaxis()->SetBinLabel(3, "met_filters_pass");
    h_cutflow->GetXaxis()->SetBinLabel(4, "pt_pass");
    h_cutflow->GetXaxis()->SetBinLabel(5, "pixel_seed_veto_pass");
    h_cutflow->GetXaxis()->SetBinLabel(6, "beamhalo_pass");
    h_cutflow->GetXaxis()->SetBinLabel(7, "tight_id_pass");
    h_cutflow->GetXaxis()->SetBinLabel(8, "met_cuts_pass");
    h_cutflow->GetXaxis()->SetBinLabel(9, "jet_veto_pass");
    h_cutflow->GetXaxis()->SetBinLabel(10," barrel_pass");
    h_cutflow->GetXaxis()->SetBinLabel(11,"barrel_MIP_pass");
    h_cutflow->GetXaxis()->SetBinLabel(12,"barrel_passSigmaIetaIeta");
    h_cutflow->GetXaxis()->SetBinLabel(13,",barrel_passSigmaIphiIphi");
    h_cutflow->GetXaxis()->SetBinLabel(14,"barrel_passEtaWidth");
    h_cutflow->GetXaxis()->SetBinLabel(15,",barrel_passTiming");
    h_cutflow->GetXaxis()->SetBinLabel(16," barrel_all_cuts");
    h_cutflow->GetXaxis()->SetBinLabel(17, "endcap_pass");

    h_pT_pho_EB->Write();
    h_pT_pho_EE->Write();
    h_eta_pho_EB->Write();
    h_eta_pho_EE->Write();
    h_phi_pho_EB->Write();
    h_phi_pho_EE->Write();
    h_pt_Met_EB->Write();
    h_phi_Met_EB->Write();
    h_pt_Met_EE->Write();
    h_phi_Met_EE->Write();
    h_jet_deltaPhi_MET->Write();
    h_deltaR_pho_jet->Write();
    h_pho_etawidth_EB->Write();
    h_pho_etawidth_EE->Write();
    h_cutflow->Write();
    
    // Force write and close
    outputFile->Write();
    outputFile->Close();
}


         
  void processMultipleFiles(const std::vector<std::string>& inputFiles, const std::string& outputFileName) {
   
    setupOutputHistograms(outputFileName);
    for (size_t i = 0; i < inputFiles.size(); i++) {      
      TFile* inputFile = TFile::Open(inputFiles[i].c_str(), "READ");
      TTree* tree = (TTree*)inputFile->Get("EventTree");
      setupInputBranches(tree);
      processEvents(tree);
      inputFile->Close();
      delete inputFile;
    }
    
    saveHistograms();
  }

  void processFileList(const std::string& fileListPath, const std::string& outputFileName) {
    std::ifstream fileList(fileListPath);
    std::vector<std::string> inputFiles;
    std::string line;
    
    while (std::getline(fileList, line)) {
      line.erase(0, line.find_first_not_of(" \t\r\n"));
      line.erase(line.find_last_not_of(" \t\r\n") + 1);
      
      if (!line.empty() && line[0] != '#') {
        inputFiles.push_back(line);
      }
    }
    fileList.close();
    processMultipleFiles(inputFiles, outputFileName);
  }
};


int main(int argc, char* argv[]) {  
  std::string outputFileName = argv[1];
  std::string fileListPath = argv[3];
  
  PhotonAnalyser analyser;
  analyser.processFileList(fileListPath, outputFileName);
   
  return 0;
}
