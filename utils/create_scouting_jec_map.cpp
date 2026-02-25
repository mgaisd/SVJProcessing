/**
 * create_scouting_jec_map.cpp
 * 
 * Fast C++ implementation for creating scouting JEC correction maps.
 * Uses ROOT TTree for I/O and explicit event loops for unique jet matching.
 * 
 * Performance: ~50-100x faster than Python implementation
 * 
 * Outputs:
 *   - ROOT file with TH2D histograms (sum, sum_squared, counts) per bin
 *   - JSON file with correction factors and uncertainties
 * 
 * Author: SVJ Scouting Analysis
 * Date: February 2026
 */

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2D.h"
#include "TError.h"
#include "TObjArray.h"

// JSON output helper
#include <sstream>

// Jet matching result
struct JetMatch {
    float scout_pt;
    float scout_eta;
    float offline_pt;
    float deltaR;
};

// Configuration
struct Config {
    std::vector<std::string> input_files;
    std::string output_dir;
    std::string jet_type;  // "Jet" or "FatJet"
    float dr_threshold;
    float pt_min;
    int max_jets;  // -1 for all jets
    
    // Binning
    int n_pt_bins;
    float pt_min_bin, pt_max_bin;
    int n_eta_bins;
    float eta_min_bin, eta_max_bin;
};

// Calculate deltaR
inline float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = phi1 - phi2;
    // Wrap phi to [-pi, pi]
    while (dphi > M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    return std::sqrt(deta*deta + dphi*dphi);
}

// Create directory (recursive)
void createDirectory(const std::string& path) {
    // Use system mkdir -p for simplicity
    std::string cmd = "mkdir -p " + path;
    int ret = system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "Warning: Could not create directory " << path << std::endl;
    }
}

// Match jets with unique matching (each offline jet matched to at most one scout)
std::vector<JetMatch> matchJets(
    const std::vector<float>& scout_pt,
    const std::vector<float>& scout_eta,
    const std::vector<float>& scout_phi,
    const std::vector<float>& offline_pt,
    const std::vector<float>& offline_eta,
    const std::vector<float>& offline_phi,
    float dr_threshold,
    int& n_scouts,
    int& n_matched)
{
    std::vector<JetMatch> matches;
    
    int n_scout = scout_pt.size();
    int n_offline = offline_pt.size();
    
    n_scouts += n_scout;
    
    if (n_scout == 0 || n_offline == 0) {
        return matches;
    }
    
    // Store all candidate matches: (scout_idx, offline_idx, deltaR)
    struct Candidate {
        int scout_idx;
        int offline_idx;
        float dr;
    };
    std::vector<Candidate> candidates;
    
    // Find best match for each scout jet
    for (int i_scout = 0; i_scout < n_scout; ++i_scout) {
        float best_dr = dr_threshold;
        int best_offline_idx = -1;
        
        for (int i_offline = 0; i_offline < n_offline; ++i_offline) {
            float dr = deltaR(scout_eta[i_scout], scout_phi[i_scout],
                            offline_eta[i_offline], offline_phi[i_offline]);
            
            if (dr < best_dr) {
                best_dr = dr;
                best_offline_idx = i_offline;
            }
        }
        
        if (best_offline_idx >= 0) {
            candidates.push_back({i_scout, best_offline_idx, best_dr});
        }
    }
    
    if (candidates.empty()) {
        return matches;
    }
    
    // Resolve conflicts: sort by offline index, then by deltaR
    std::sort(candidates.begin(), candidates.end(),
        [](const Candidate& a, const Candidate& b) {
            if (a.offline_idx != b.offline_idx) return a.offline_idx < b.offline_idx;
            return a.dr < b.dr;
        });
    
    // Keep only first occurrence of each unique offline jet (smallest deltaR)
    std::vector<bool> offline_used(n_offline, false);
    
    for (const auto& cand : candidates) {
        if (!offline_used[cand.offline_idx]) {
            offline_used[cand.offline_idx] = true;
            
            // Sanity checks
            float s_pt = scout_pt[cand.scout_idx];
            float o_pt = offline_pt[cand.offline_idx];
            
            if (s_pt > 0 && o_pt > 0 && std::isfinite(s_pt) && std::isfinite(o_pt)) {
                JetMatch match;
                match.scout_pt = s_pt;
                match.scout_eta = scout_eta[cand.scout_idx];
                match.offline_pt = o_pt;
                match.deltaR = cand.dr;
                matches.push_back(match);
                n_matched++;
            }
        }
    }
    
    return matches;
}

// Process one ROOT file for BOTH jet types simultaneously
void processFileBoth(
    const std::string& filename,
    const Config& config_jet, TH2D* h_sum_jet, TH2D* h_sum_sq_jet, TH2D* h_counts_jet,
    long& n_total_jets, long& n_rejected_neg_jet, long& n_rejected_outliers_jet,
    const Config& config_fatjet, TH2D* h_sum_fatjet, TH2D* h_sum_sq_fatjet, TH2D* h_counts_fatjet,
    long& n_total_fatjets, long& n_rejected_neg_fatjet, long& n_rejected_outliers_fatjet)
{
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "  Warning: Cannot open file: " << filename << std::endl;
        return;
    }
    
    TTree* tree = (TTree*)file->Get("mmtree/Events");
    if (!tree) {
        std::cerr << "  Warning: Cannot find tree mmtree/Events in " << filename << std::endl;
        file->Close();
        return;
    }
    
    // Enable TTreeCache for faster reading (30 MB cache)
    tree->SetCacheSize(30000000);
    tree->AddBranchToCache("*", kTRUE);
    tree->StopCacheLearningPhase();
    
    // Set up branches for BOTH jet types
    std::vector<float> *jet_pt = nullptr, *jet_eta = nullptr, *jet_phi = nullptr;
    std::vector<float> *offline_jet_pt = nullptr, *offline_jet_eta = nullptr, *offline_jet_phi = nullptr;
    std::vector<float> *fatjet_pt = nullptr, *fatjet_eta = nullptr, *fatjet_phi = nullptr;
    std::vector<float> *offline_fatjet_pt = nullptr, *offline_fatjet_eta = nullptr, *offline_fatjet_phi = nullptr;
    
    tree->SetBranchAddress("Jet_pt", &jet_pt);
    tree->SetBranchAddress("Jet_eta", &jet_eta);
    tree->SetBranchAddress("Jet_phi", &jet_phi);
    tree->SetBranchAddress("OfflineJet_pt", &offline_jet_pt);
    tree->SetBranchAddress("OfflineJet_eta", &offline_jet_eta);
    tree->SetBranchAddress("OfflineJet_phi", &offline_jet_phi);
    
    tree->SetBranchAddress("FatJet_pt", &fatjet_pt);
    tree->SetBranchAddress("FatJet_eta", &fatjet_eta);
    tree->SetBranchAddress("FatJet_phi", &fatjet_phi);
    tree->SetBranchAddress("OfflineFatJet_pt", &offline_fatjet_pt);
    tree->SetBranchAddress("OfflineFatJet_eta", &offline_fatjet_eta);
    tree->SetBranchAddress("OfflineFatJet_phi", &offline_fatjet_phi);
    
    int n_events = tree->GetEntries();
    
    // Process events - fill BOTH jet types simultaneously
    for (int i_evt = 0; i_evt < n_events; ++i_evt) {
        tree->GetEntry(i_evt);
        
        // ===== Process AK4 (Jet) =====
        if (jet_pt && offline_jet_pt) {
            std::vector<float> s_pt, s_eta, s_phi;
            std::vector<float> o_pt, o_eta, o_phi;
            
            int max_scout = (config_jet.max_jets > 0) ? std::min(config_jet.max_jets, (int)jet_pt->size()) : jet_pt->size();
            
            for (int i = 0; i < max_scout; ++i) {
                if ((*jet_pt)[i] > config_jet.pt_min) {
                    s_pt.push_back((*jet_pt)[i]);
                    s_eta.push_back((*jet_eta)[i]);
                    s_phi.push_back((*jet_phi)[i]);
                }
            }
            
            for (size_t i = 0; i < offline_jet_pt->size(); ++i) {
                if ((*offline_jet_pt)[i] > config_jet.pt_min) {
                    o_pt.push_back((*offline_jet_pt)[i]);
                    o_eta.push_back((*offline_jet_eta)[i]);
                    o_phi.push_back((*offline_jet_phi)[i]);
                }
            }
            
            int n_scouts_temp = 0, n_matched_temp = 0;
            std::vector<JetMatch> matches = matchJets(s_pt, s_eta, s_phi, o_pt, o_eta, o_phi,
                                                      config_jet.dr_threshold, n_scouts_temp, n_matched_temp);
            
            for (const auto& match : matches) {
                n_total_jets++;
                float ratio = match.offline_pt / match.scout_pt;
                
                if (ratio <= 0.1 || ratio >= 10.0) {
                    n_rejected_outliers_jet++;
                    continue;
                }
                
                int bin_x = h_sum_jet->GetXaxis()->FindBin(match.scout_pt);
                int bin_y = h_sum_jet->GetYaxis()->FindBin(match.scout_eta);
                
                if (bin_x >= 1 && bin_x <= h_sum_jet->GetNbinsX() &&
                    bin_y >= 1 && bin_y <= h_sum_jet->GetNbinsY()) {
                    h_sum_jet->Fill(match.scout_pt, match.scout_eta, ratio);
                    h_sum_sq_jet->Fill(match.scout_pt, match.scout_eta, ratio * ratio);
                    h_counts_jet->Fill(match.scout_pt, match.scout_eta);
                }
            }
        }
        
        // ===== Process AK8 (FatJet) =====
        if (fatjet_pt && offline_fatjet_pt) {
            std::vector<float> s_pt, s_eta, s_phi;
            std::vector<float> o_pt, o_eta, o_phi;
            
            int max_scout = (config_fatjet.max_jets > 0) ? std::min(config_fatjet.max_jets, (int)fatjet_pt->size()) : fatjet_pt->size();
            
            for (int i = 0; i < max_scout; ++i) {
                if ((*fatjet_pt)[i] > config_fatjet.pt_min) {
                    s_pt.push_back((*fatjet_pt)[i]);
                    s_eta.push_back((*fatjet_eta)[i]);
                    s_phi.push_back((*fatjet_phi)[i]);
                }
            }
            
            for (size_t i = 0; i < offline_fatjet_pt->size(); ++i) {
                if ((*offline_fatjet_pt)[i] > config_fatjet.pt_min) {
                    o_pt.push_back((*offline_fatjet_pt)[i]);
                    o_eta.push_back((*offline_fatjet_eta)[i]);
                    o_phi.push_back((*offline_fatjet_phi)[i]);
                }
            }
            
            int n_scouts_temp = 0, n_matched_temp = 0;
            std::vector<JetMatch> matches = matchJets(s_pt, s_eta, s_phi, o_pt, o_eta, o_phi,
                                                      config_fatjet.dr_threshold, n_scouts_temp, n_matched_temp);
            
            for (const auto& match : matches) {
                n_total_fatjets++;
                float ratio = match.offline_pt / match.scout_pt;
                
                if (ratio <= 0.1 || ratio >= 10.0) {
                    n_rejected_outliers_fatjet++;
                    continue;
                }
                
                int bin_x = h_sum_fatjet->GetXaxis()->FindBin(match.scout_pt);
                int bin_y = h_sum_fatjet->GetYaxis()->FindBin(match.scout_eta);
                
                if (bin_x >= 1 && bin_x <= h_sum_fatjet->GetNbinsX() &&
                    bin_y >= 1 && bin_y <= h_sum_fatjet->GetNbinsY()) {
                    h_sum_fatjet->Fill(match.scout_pt, match.scout_eta, ratio);
                    h_sum_sq_fatjet->Fill(match.scout_pt, match.scout_eta, ratio * ratio);
                    h_counts_fatjet->Fill(match.scout_pt, match.scout_eta);
                }
            }
        }
    }
    
    file->Close();
}

// Save results to JSON
void saveJSON(const Config& config, TH2D* h_sum, TH2D* h_sum_sq, TH2D* h_counts,
              long n_total, long n_rejected_neg, long n_rejected_outliers)
{
    std::string json_file = config.output_dir + "/scouting_jec_corrections_" + config.jet_type + ".json";
    std::ofstream out(json_file);
    
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create JSON file: " << json_file << std::endl;
        return;
    }
    
    out << std::fixed << std::setprecision(6);
    out << "{\n";
    out << "  \"jet_type\": \"" << config.jet_type << "\",\n";
    out << "  \"dr_threshold\": " << config.dr_threshold << ",\n";
    out << "  \"pt_min\": " << config.pt_min << ",\n";
    
    // Write bin edges
    out << "  \"pt_bins\": [";
    for (int i = 1; i <= h_sum->GetNbinsX() + 1; ++i) {
        out << h_sum->GetXaxis()->GetBinLowEdge(i);
        if (i <= h_sum->GetNbinsX()) out << ", ";
    }
    out << "],\n";
    
    out << "  \"eta_bins\": [";
    for (int i = 1; i <= h_sum->GetNbinsY() + 1; ++i) {
        out << h_sum->GetYaxis()->GetBinLowEdge(i);
        if (i <= h_sum->GetNbinsY()) out << ", ";
    }
    out << "],\n";
    
    // Calculate corrections and uncertainties
    const int min_count = 10;
    out << "  \"corrections\": [\n";
    for (int ix = 1; ix <= h_sum->GetNbinsX(); ++ix) {
        out << "    [";
        for (int iy = 1; iy <= h_sum->GetNbinsY(); ++iy) {
            double sum = h_sum->GetBinContent(ix, iy);
            double count = h_counts->GetBinContent(ix, iy);
            double corr = 1.0;
            if (count >= min_count) {
                corr = sum / count;
                if (corr <= 0.1 || corr >= 10.0) corr = 1.0;
            }
            out << corr;
            if (iy < h_sum->GetNbinsY()) out << ", ";
        }
        out << "]";
        if (ix < h_sum->GetNbinsX()) out << ",";
        out << "\n";
    }
    out << "  ],\n";

    out << "  \"uncertainties\": [\n";
    for (int ix = 1; ix <= h_sum->GetNbinsX(); ++ix) {
        out << "    [";
        for (int iy = 1; iy <= h_sum->GetNbinsY(); ++iy) {
            double sum = h_sum->GetBinContent(ix, iy);
            double sum_sq = h_sum_sq->GetBinContent(ix, iy);
            double count = h_counts->GetBinContent(ix, iy);
            double unc = 0.0;
            if (count >= min_count) {
                if (count > 1) {
                    double mean = sum / count;
                    double mean_sq = sum_sq / count;
                    double variance = mean_sq - mean * mean;
                    unc = (variance > 0) ? std::sqrt(variance) : 0.0;
                }
            }
            out << unc;
            if (iy < h_sum->GetNbinsY()) out << ", ";
        }
        out << "]";
        if (ix < h_sum->GetNbinsX()) out << ",";
        out << "\n";
    }
    out << "  ],\n";

    out << "  \"stat_uncertainties\": [\n";
    for (int ix = 1; ix <= h_sum->GetNbinsX(); ++ix) {
        out << "    [";
        for (int iy = 1; iy <= h_sum->GetNbinsY(); ++iy) {
            double sum = h_sum->GetBinContent(ix, iy);
            double sum_sq = h_sum_sq->GetBinContent(ix, iy);
            double count = h_counts->GetBinContent(ix, iy);
            double stat_unc = 0.0;
            if (count >= min_count) {
                if (count > 1) {
                    double mean = sum / count;
                    double mean_sq = sum_sq / count;
                    double variance = mean_sq - mean * mean;
                    double unc = (variance > 0) ? std::sqrt(variance) : 0.0;
                    stat_unc = unc / std::sqrt(count);
                }
            }
            out << stat_unc;
            if (iy < h_sum->GetNbinsY()) out << ", ";
        }
        out << "]";
        if (ix < h_sum->GetNbinsX()) out << ",";
        out << "\n";
    }
    out << "  ],\n";
    
    out << "  \"counts\": [\n";
    for (int ix = 1; ix <= h_sum->GetNbinsX(); ++ix) {
        out << "    [";
        for (int iy = 1; iy <= h_sum->GetNbinsY(); ++iy) {
            out << (int)h_counts->GetBinContent(ix, iy);
            if (iy < h_sum->GetNbinsY()) out << ", ";
        }
        out << "]";
        if (ix < h_sum->GetNbinsX()) out << ",";
        out << "\n";
    }
    out << "  ],\n";
    
    // Statistics
    long n_accepted = n_total - n_rejected_neg - n_rejected_outliers;
    out << "  \"statistics\": {\n";
    out << "    \"total_matched\": " << n_total << ",\n";
    out << "    \"rejected_negative_pt\": " << n_rejected_neg << ",\n";
    out << "    \"rejected_outliers\": " << n_rejected_outliers << ",\n";
    out << "    \"accepted\": " << n_accepted << "\n";
    out << "  }\n";
    out << "}\n";
    
    out.close();
    std::cout << "Saved JSON: " << json_file << std::endl;
}

int main(int argc, char** argv) {
    // Suppress ROOT warnings about TClassTable
    gErrorIgnoreLevel = kError;
    
    std::cout << "=============================================================\n";
    std::cout << "Scouting JEC Map Creator (C++ version)\n";
    std::cout << "=============================================================\n";
    
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <file1.root> [file2.root ...]\n";
        std::cerr << "  Processes both Jet (AK4) and FatJet (AK8) in a single pass\n";
        return 1;
    }
    
    // Get input files
    std::vector<std::string> input_files;
    for (int i = 1; i < argc; ++i) {
        input_files.push_back(argv[i]);
    }
    
    // Configuration for AK4 (Jet)
    Config config_jet;
    config_jet.jet_type = "Jet";
    config_jet.dr_threshold = 0.15;
    config_jet.output_dir = "test";
    config_jet.pt_min = 15.0;
    config_jet.pt_min_bin = 15.0;
    config_jet.pt_max_bin = 2000.0;
    config_jet.n_pt_bins = 30;
    config_jet.eta_min_bin = -5.0;
    config_jet.eta_max_bin = 5.0;
    config_jet.n_eta_bins = 24;
    config_jet.max_jets = -1;  // All jets
    config_jet.input_files = input_files;
    
    // Configuration for AK8 (FatJet)
    Config config_fatjet;
    config_fatjet.jet_type = "FatJet";
    config_fatjet.dr_threshold = 0.15;
    config_fatjet.output_dir = "test";
    config_fatjet.pt_min = 150.0;
    config_fatjet.pt_min_bin = 150.0;
    config_fatjet.pt_max_bin = 2000.0;
    config_fatjet.n_pt_bins = 30;
    config_fatjet.eta_min_bin = -2.4;
    config_fatjet.eta_max_bin = 2.4;
    config_fatjet.n_eta_bins = 24;
    config_fatjet.max_jets = 2;  // Leading 2 jets
    config_fatjet.input_files = input_files;
    
    std::cout << "Configuration:\n";
    std::cout << "  Processing: Jet (AK4) + FatJet (AK8) simultaneously\n";
    std::cout << "  deltaR threshold: " << config_jet.dr_threshold << "\n";
    std::cout << "  Input files: " << input_files.size() << "\n";
    std::cout << "  Output directory: " << config_jet.output_dir << "\n\n";
    
    // Create output directory
    createDirectory(config_jet.output_dir);
    
    // Create histograms for AK4 (Jet)
    TH2D* h_sum_jet = new TH2D("h_sum_jet", "Sum of ratios (Jet)",
                          config_jet.n_pt_bins, config_jet.pt_min_bin, config_jet.pt_max_bin,
                          config_jet.n_eta_bins, config_jet.eta_min_bin, config_jet.eta_max_bin);
    TH2D* h_sum_sq_jet = new TH2D("h_sum_sq_jet", "Sum of ratios squared (Jet)",
                             config_jet.n_pt_bins, config_jet.pt_min_bin, config_jet.pt_max_bin,
                             config_jet.n_eta_bins, config_jet.eta_min_bin, config_jet.eta_max_bin);
    TH2D* h_counts_jet = new TH2D("h_counts_jet", "Number of jets (Jet)",
                             config_jet.n_pt_bins, config_jet.pt_min_bin, config_jet.pt_max_bin,
                             config_jet.n_eta_bins, config_jet.eta_min_bin, config_jet.eta_max_bin);
    
    h_sum_jet->Sumw2();
    h_sum_sq_jet->Sumw2();
    h_counts_jet->Sumw2();
    
    // Create histograms for AK8 (FatJet)
    TH2D* h_sum_fatjet = new TH2D("h_sum_fatjet", "Sum of ratios (FatJet)",
                          config_fatjet.n_pt_bins, config_fatjet.pt_min_bin, config_fatjet.pt_max_bin,
                          config_fatjet.n_eta_bins, config_fatjet.eta_min_bin, config_fatjet.eta_max_bin);
    TH2D* h_sum_sq_fatjet = new TH2D("h_sum_sq_fatjet", "Sum of ratios squared (FatJet)",
                             config_fatjet.n_pt_bins, config_fatjet.pt_min_bin, config_fatjet.pt_max_bin,
                             config_fatjet.n_eta_bins, config_fatjet.eta_min_bin, config_fatjet.eta_max_bin);
    TH2D* h_counts_fatjet = new TH2D("h_counts_fatjet", "Number of jets (FatJet)",
                             config_fatjet.n_pt_bins, config_fatjet.pt_min_bin, config_fatjet.pt_max_bin,
                             config_fatjet.n_eta_bins, config_fatjet.eta_min_bin, config_fatjet.eta_max_bin);
    
    h_sum_fatjet->Sumw2();
    h_sum_sq_fatjet->Sumw2();
    h_counts_fatjet->Sumw2();
    
    // Statistics for both jet types
    long n_total_jets = 0, n_rejected_neg = 0, n_rejected_outliers = 0;
    long n_total_fatjets = 0, n_rejected_neg_fatjet = 0, n_rejected_outliers_fatjet = 0;
    
    // Process files - BOTH jet types simultaneously
    std::cout << "Processing files (both Jet and FatJet)...\n";
    for (size_t i = 0; i < input_files.size(); ++i) {
        double pct = 100.0 * (i + 1) / input_files.size();
        std::cout << "  File " << (i+1) << "/" << input_files.size() 
                  << " (" << std::fixed << std::setprecision(1) << pct << "%)\n";
        
        // Process BOTH jet types from same file in one pass
        processFileBoth(input_files[i],
                       config_jet, h_sum_jet, h_sum_sq_jet, h_counts_jet,
                       n_total_jets, n_rejected_neg, n_rejected_outliers,
                       config_fatjet, h_sum_fatjet, h_sum_sq_fatjet, h_counts_fatjet,
                       n_total_fatjets, n_rejected_neg_fatjet, n_rejected_outliers_fatjet);
    }
    
    // Print statistics for AK4 (Jet)
    std::cout << "\n============================================\n";
    std::cout << "AK4 (Jet) Summary:\n";
    std::cout << "  Total matched jets: " << n_total_jets << "\n";
    std::cout << "  Rejected (outliers): " << n_rejected_outliers 
              << " (" << (n_total_jets > 0 ? 100.0*n_rejected_outliers/n_total_jets : 0) << "%)\n";
    std::cout << "  Accepted: " << (n_total_jets - n_rejected_neg - n_rejected_outliers)
              << " (" << (n_total_jets > 0 ? 100.0*(n_total_jets-n_rejected_neg-n_rejected_outliers)/n_total_jets : 0) << "%)\n";


    // Print statistics for AK8 (FatJet)
    std::cout << "\nAK8 (FatJet) Summary:\n";
    std::cout << "  Total matched jets: " << n_total_fatjets << "\n";
    std::cout << "  Rejected (outliers): " << n_rejected_outliers_fatjet 
              << " (" << (n_total_fatjets > 0 ? 100.0*n_rejected_outliers_fatjet/n_total_fatjets : 0) << "%)\n";
    std::cout << "  Accepted: " << (n_total_fatjets - n_rejected_neg_fatjet - n_rejected_outliers_fatjet)
              << " (" << (n_total_fatjets > 0 ? 100.0*(n_total_fatjets-n_rejected_neg_fatjet-n_rejected_outliers_fatjet)/n_total_fatjets : 0) << "%)\n";

    
    // Save ROOT histograms for AK4
    std::string root_file_jet = config_jet.output_dir + "/scouting_jec_histograms_Jet.root";
    TFile* out_file_jet = new TFile(root_file_jet.c_str(), "RECREATE");
    h_sum_jet->Write("h_sum");
    h_sum_sq_jet->Write("h_sum_sq");
    h_counts_jet->Write("h_counts");
    out_file_jet->Close();
    std::cout << "\nSaved ROOT file: " << root_file_jet << std::endl;
    
    // Save ROOT histograms for AK8
    std::string root_file_fatjet = config_fatjet.output_dir + "/scouting_jec_histograms_FatJet.root";
    TFile* out_file_fatjet = new TFile(root_file_fatjet.c_str(), "RECREATE");
    h_sum_fatjet->Write("h_sum");
    h_sum_sq_fatjet->Write("h_sum_sq");
    h_counts_fatjet->Write("h_counts");
    out_file_fatjet->Close();
    std::cout << "Saved ROOT file: " << root_file_fatjet << std::endl;
    
    // Save JSON for AK4
    saveJSON(config_jet, h_sum_jet, h_sum_sq_jet, h_counts_jet,
             n_total_jets, n_rejected_neg, n_rejected_outliers);
    
    // Save JSON for AK8
    saveJSON(config_fatjet, h_sum_fatjet, h_sum_sq_fatjet, h_counts_fatjet,
             n_total_fatjets, n_rejected_neg_fatjet, n_rejected_outliers_fatjet);
    
    // Cleanup
    delete h_sum_jet;
    delete h_sum_sq_jet;
    delete h_counts_jet;
    delete h_sum_fatjet;
    delete h_sum_sq_fatjet;
    delete h_counts_fatjet;
    
    std::cout << "\n=============================================================\n";
    std::cout << "Done! Run plot_jec_corrections.py to create plots.\n";
    std::cout << "=============================================================\n";
    
    return 0;
}
