#!/usr/bin/env python3
"""
Diagnostic script to check JEC corrections on FatJets
Compares gen, corrected (nominal), and uncorrected jet pT distributions
Also plots response (reco/gen) in bins of gen pT
"""

import awkward as ak
import uproot
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from pathlib import Path

# QCD sample configuration
# Adjust paths and cross-sections as needed
QCD_SAMPLES = {
    'QCD_HT300to500': {
        'path': '/ceph/mgais/Run2ScoutingSkims_JEC/2018/QCD_HT300to500/nominal/*.root',
        'xsec': 323400.0,  # pb
        'color': 'tab:blue'
    },
    'QCD_HT500to700': {
        'path': '/ceph/mgais/Run2ScoutingSkims_JEC/2018/QCD_HT500to700/nominal/*.root',
        'xsec': 30140.0,
        'color': 'tab:orange'
    },
    'QCD_HT700to1000': {
        'path': '/ceph/mgais/Run2ScoutingSkims_JEC/2018/QCD_HT700to1000/nominal/*.root',
        'xsec': 6310.0,
        'color': 'tab:green'
    },
    'QCD_HT1000to1500': {
        'path': '/ceph/mgais/Run2ScoutingSkims_JEC/2018/QCD_HT1000to1500/nominal/*.root',
        'xsec': 1094.0,
        'color': 'tab:red'
    },
    'QCD_HT1500to2000': {
        'path': '/ceph/mgais/Run2ScoutingSkims_JEC/2018/QCD_HT1500to2000/nominal/*.root',
        'xsec': 99.38,
        'color': 'tab:purple'
    },
    'QCD_HT2000toInf': {
        'path': '/ceph/mgais/Run2ScoutingSkims_JEC/2018/QCD_HT2000toInf/nominal/*.root',
        'xsec': 20.20,
        'color': 'tab:brown'
    },
}
LUMI = 59830.0  # pb^-1 for 2018


def load_sample(sample_name, sample_info, max_files=1):
    """Load events from a sample"""
    print(f"\nLoading {sample_name}...")
    
    file_pattern = sample_info['path']
    files = sorted(glob.glob(file_pattern))
    
    if not files:
        print(f"  WARNING: No files found for {sample_name} at {file_pattern}")
        return None
    
    print(f"  Found {len(files)} files, loading first {max_files}")
    files = files[:max_files]
    
    events_list = []
    total_gen_weight_sum = 0
    
    for fpath in files:
        try:
            with uproot.open(fpath) as f:
                tree = f['Events']
                
                # Load branches we need
                branches = [
                    'FatJet_pt', 'FatJet_eta', 'FatJet_phi', 'FatJet_mass',
                    'FatJet_pt_uncorr', 'FatJet_mass_uncorr',
                    'GenFatJet_pt', 'GenFatJet_eta', 'GenFatJet_phi', 'GenFatJet_mass',
                    'genWeight'
                ]
                
                # Check if JEC variations exist
                all_branches = tree.keys()
                has_jec_up = 'FatJet_pt_jec_up' in all_branches
                has_jec_down = 'FatJet_pt_jec_down' in all_branches
                
                if has_jec_up:
                    branches.extend(['FatJet_pt_jec_up', 'FatJet_mass_jec_up'])
                if has_jec_down:
                    branches.extend(['FatJet_pt_jec_down', 'FatJet_mass_jec_down'])
                
                events = tree.arrays(branches, library='ak')
                events_list.append(events)
                
                # Sum genWeights for normalization
                total_gen_weight_sum += ak.sum(events.genWeight)
                
        except Exception as e:
            print(f"  Error loading {fpath}: {e}")
            continue
    
    if not events_list:
        return None
    
    events = ak.concatenate(events_list, axis=0)
    
    # Calculate event weights
    xsec = sample_info['xsec']
    weight = xsec * LUMI / total_gen_weight_sum
    events['weight'] = events.genWeight * weight
    
    # Store metadata
    events.has_jec_up = has_jec_up
    events.has_jec_down = has_jec_down
    
    print(f"  Loaded {len(events)} events")
    print(f"  XSec weight factor: {weight:.6f}")
    print(f"  JEC variations: up={has_jec_up}, down={has_jec_down}")
    
    return events


def match_jets(reco_jets, gen_jets, dr_max=0.4):
    """Match reco jets to gen jets using dR matching"""
    # Broadcast for pairwise distance calculation
    reco_eta = ak.unflatten(reco_jets.eta, ak.num(reco_jets.eta))
    reco_phi = ak.unflatten(reco_jets.phi, ak.num(reco_jets.phi))
    gen_eta = ak.unflatten(gen_jets.eta, ak.num(gen_jets.eta))
    gen_phi = ak.unflatten(gen_jets.phi, ak.num(gen_jets.phi))
    
    # Calculate dR between all reco-gen pairs
    deta = reco_eta[:, :, np.newaxis] - gen_eta[:, np.newaxis, :]
    dphi = np.abs(reco_phi[:, :, np.newaxis] - gen_phi[:, np.newaxis, :])
    dphi = np.where(dphi > np.pi, 2*np.pi - dphi, dphi)
    dr = np.sqrt(deta**2 + dphi**2)
    
    # Find closest gen jet for each reco jet
    min_dr_idx = ak.argmin(dr, axis=2, keepdims=False)
    min_dr = ak.min(dr, axis=2)
    
    # Mask for matched jets (dR < dr_max)
    matched_mask = min_dr < dr_max
    
    return matched_mask, min_dr_idx


def plot_pt_distributions(all_events, output_dir='plots'):
    """Plot pT distributions for gen, corrected, and uncorrected jets"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
    
    bins = np.linspace(0, 2000, 51)
    
    # Collect data from all samples
    gen_pt_all = []
    corr_pt_all = []
    uncorr_pt_all = []
    weights_gen = []
    weights_corr = []
    weights_uncorr = []
    
    for sample_name, events in all_events.items():
        if events is None:
            continue
            
        weight = events.weight
        
        # Gen jets
        gen_pt = ak.flatten(events.GenFatJet_pt)
        gen_weight = ak.broadcast_arrays(weight, events.GenFatJet_pt)[0]
        gen_weight = ak.flatten(gen_weight)
        
        gen_pt_all.append(ak.to_numpy(gen_pt))
        weights_gen.append(ak.to_numpy(gen_weight))
        
        # Corrected jets (nominal)
        corr_pt = ak.flatten(events.FatJet_pt)
        corr_weight = ak.broadcast_arrays(weight, events.FatJet_pt)[0]
        corr_weight = ak.flatten(corr_weight)
        
        corr_pt_all.append(ak.to_numpy(corr_pt))
        weights_corr.append(ak.to_numpy(corr_weight))
        
        # Uncorrected jets
        uncorr_pt = ak.flatten(events.FatJet_pt_raw)
        uncorr_weight = ak.broadcast_arrays(weight, events.FatJet_pt_raw)[0]
        uncorr_weight = ak.flatten(uncorr_weight)
        
        uncorr_pt_all.append(ak.to_numpy(uncorr_pt))
        weights_uncorr.append(ak.to_numpy(uncorr_weight))
    
    # Concatenate all samples
    gen_pt_all = np.concatenate(gen_pt_all)
    corr_pt_all = np.concatenate(corr_pt_all)
    uncorr_pt_all = np.concatenate(uncorr_pt_all)
    weights_gen = np.concatenate(weights_gen)
    weights_corr = np.concatenate(weights_corr)
    weights_uncorr = np.concatenate(weights_uncorr)
    
    # Plot distributions
    ax1.hist(gen_pt_all, bins=bins, weights=weights_gen, 
             histtype='step', linewidth=2, label='Gen', color='black')
    ax1.hist(corr_pt_all, bins=bins, weights=weights_corr,
             histtype='step', linewidth=2, label='Corrected (nominal)', color='tab:blue')
    ax1.hist(uncorr_pt_all, bins=bins, weights=weights_uncorr,
             histtype='step', linewidth=2, label='Uncorrected (raw)', color='tab:red', linestyle='--')
    
    ax1.set_ylabel('Events', fontsize=12)
    ax1.legend(fontsize=11)
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)
    ax1.text(0.05, 0.95, 'AK8 Jets', transform=ax1.transAxes, 
             fontsize=11, verticalalignment='top', 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Ratio plot
    gen_hist, _ = np.histogram(gen_pt_all, bins=bins, weights=weights_gen)
    corr_hist, _ = np.histogram(corr_pt_all, bins=bins, weights=weights_corr)
    uncorr_hist, _ = np.histogram(uncorr_pt_all, bins=bins, weights=weights_uncorr)
    
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio_corr = np.where(gen_hist > 0, corr_hist / gen_hist, np.nan)
        ratio_uncorr = np.where(gen_hist > 0, uncorr_hist / gen_hist, np.nan)
    
    ax2.plot(bin_centers, ratio_corr, 'o-', color='tab:blue', label='Corrected / Gen')
    ax2.plot(bin_centers, ratio_uncorr, 's--', color='tab:red', label='Uncorrected / Gen')
    ax2.axhline(1.0, color='black', linestyle=':', alpha=0.5)
    
    ax2.set_xlabel('Jet $p_T$ [GeV]', fontsize=12)
    ax2.set_ylabel('Ratio to Gen', fontsize=12)
    ax2.legend(fontsize=11)
    ax2.set_ylim(0.5, 1.5)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fatjet_pt_comparison.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{output_dir}/fatjet_pt_comparison.pdf', bbox_inches='tight')
    print(f"\nSaved pT comparison plots to {output_dir}/fatjet_pt_comparison.{{png,pdf}}")
    plt.close()


def plot_response(all_events, output_dir='plots'):
    """Plot response (reco/gen) vs gen pT for corrected and uncorrected jets"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Define gen pT bins
    gen_pt_bins = np.array([0, 300, 400, 500, 600, 800, 1000, 1500, 2000])
    
    # Storage for response calculations
    responses_corr = {i: [] for i in range(len(gen_pt_bins)-1)}
    responses_uncorr = {i: [] for i in range(len(gen_pt_bins)-1)}
    weights_corr = {i: [] for i in range(len(gen_pt_bins)-1)}
    weights_uncorr = {i: [] for i in range(len(gen_pt_bins)-1)}
    
    for sample_name, events in all_events.items():
        if events is None:
            continue
        
        print(f"\nProcessing {sample_name} for response calculation...")
        
        # Match reco to gen jets
        reco_jets = ak.zip({
            'pt': events.FatJet_pt,
            'eta': events.FatJet_eta,
            'phi': events.FatJet_phi,
            'pt_raw': events.FatJet_pt_raw
        })
        
        gen_jets = ak.zip({
            'pt': events.GenFatJet_pt,
            'eta': events.GenFatJet_eta,
            'phi': events.GenFatJet_phi
        })
        
        matched_mask, matched_idx = match_jets(reco_jets, gen_jets, dr_max=0.4)
        
        # Get matched jets
        matched_reco = reco_jets[matched_mask]
        matched_gen_idx_flat = ak.flatten(matched_idx[matched_mask])
        matched_gen = gen_jets[ak.unflatten(matched_gen_idx_flat, ak.num(matched_reco.pt))]
        
        # Calculate responses
        gen_pt = ak.flatten(matched_gen.pt)
        corr_pt = ak.flatten(matched_reco.pt)
        uncorr_pt = ak.flatten(matched_reco.pt_raw)
        
        response_corr = corr_pt / gen_pt
        response_uncorr = uncorr_pt / gen_pt
        
        # Event weights
        weight = events.weight
        event_weight = ak.broadcast_arrays(weight, matched_reco.pt)[0]
        event_weight = ak.flatten(event_weight)
        
        # Fill into pT bins
        gen_pt_np = ak.to_numpy(gen_pt)
        response_corr_np = ak.to_numpy(response_corr)
        response_uncorr_np = ak.to_numpy(response_uncorr)
        event_weight_np = ak.to_numpy(event_weight)
        
        for i in range(len(gen_pt_bins)-1):
            mask = (gen_pt_np >= gen_pt_bins[i]) & (gen_pt_np < gen_pt_bins[i+1])
            if np.sum(mask) > 0:
                responses_corr[i].append(response_corr_np[mask])
                responses_uncorr[i].append(response_uncorr_np[mask])
                weights_corr[i].append(event_weight_np[mask])
                weights_uncorr[i].append(event_weight_np[mask])
    
    # Calculate mean and std in each bin
    bin_centers = (gen_pt_bins[:-1] + gen_pt_bins[1:]) / 2
    bin_widths = (gen_pt_bins[1:] - gen_pt_bins[:-1]) / 2
    
    mean_response_corr = []
    std_response_corr = []
    mean_response_uncorr = []
    std_response_uncorr = []
    
    for i in range(len(gen_pt_bins)-1):
        if responses_corr[i]:
            resp_c = np.concatenate(responses_corr[i])
            w_c = np.concatenate(weights_corr[i])
            mean_response_corr.append(np.average(resp_c, weights=w_c))
            std_response_corr.append(np.sqrt(np.average((resp_c - mean_response_corr[-1])**2, weights=w_c)))
        else:
            mean_response_corr.append(np.nan)
            std_response_corr.append(np.nan)
            
        if responses_uncorr[i]:
            resp_u = np.concatenate(responses_uncorr[i])
            w_u = np.concatenate(weights_uncorr[i])
            mean_response_uncorr.append(np.average(resp_u, weights=w_u))
            std_response_uncorr.append(np.sqrt(np.average((resp_u - mean_response_uncorr[-1])**2, weights=w_u)))
        else:
            mean_response_uncorr.append(np.nan)
            std_response_uncorr.append(np.nan)
    
    mean_response_corr = np.array(mean_response_corr)
    std_response_corr = np.array(std_response_corr)
    mean_response_uncorr = np.array(mean_response_uncorr)
    std_response_uncorr = np.array(std_response_uncorr)
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.errorbar(bin_centers, mean_response_corr, yerr=std_response_corr, xerr=bin_widths,
                fmt='o-', color='tab:blue', label='Corrected', capsize=5, markersize=8)
    ax.errorbar(bin_centers+20, mean_response_uncorr, yerr=std_response_uncorr, xerr=bin_widths,
                fmt='s--', color='tab:red', label='Uncorrected (raw)', capsize=5, markersize=8)
    
    ax.axhline(1.0, color='black', linestyle=':', alpha=0.5, label='Perfect response')
    
    ax.set_xlabel('Gen Jet $p_T$ [GeV]', fontsize=13)
    ax.set_ylabel('Response (Reco $p_T$ / Gen $p_T$)', fontsize=13)
    ax.set_ylim(0.7, 1.15)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.text(0.05, 0.95, 'AK8 Jets\n$\Delta R <$ 0.4 matching', 
            transform=ax.transAxes, fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fatjet_response_vs_genpt.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{output_dir}/fatjet_response_vs_genpt.pdf', bbox_inches='tight')
    print(f"Saved response plots to {output_dir}/fatjet_response_vs_genpt.{{png,pdf}}")
    plt.close()
    
    # Print summary table
    print("\n" + "="*80)
    print("Response Summary (Mean ± Std)")
    print("="*80)
    print(f"{'Gen pT Bin [GeV]':<20} {'Corrected':<25} {'Uncorrected':<25}")
    print("-"*80)
    for i in range(len(gen_pt_bins)-1):
        bin_label = f"{gen_pt_bins[i]:.0f}-{gen_pt_bins[i+1]:.0f}"
        corr_str = f"{mean_response_corr[i]:.3f} ± {std_response_corr[i]:.3f}" if not np.isnan(mean_response_corr[i]) else "N/A"
        uncorr_str = f"{mean_response_uncorr[i]:.3f} ± {std_response_uncorr[i]:.3f}" if not np.isnan(mean_response_uncorr[i]) else "N/A"
        print(f"{bin_label:<20} {corr_str:<25} {uncorr_str:<25}")
    print("="*80)


def main():
    print("="*80)
    print("JEC Diagnostic Tool for FatJets")
    print("="*80)
    
    # Load samples
    all_events = {}
    for sample_name, sample_info in QCD_SAMPLES.items():
        events = load_sample(sample_name, sample_info, max_files=1)
        all_events[sample_name] = events
    
    # Check if any samples loaded
    valid_samples = [s for s, e in all_events.items() if e is not None]
    if not valid_samples:
        print("\nERROR: No samples loaded successfully!")
        return
    
    print(f"\nSuccessfully loaded {len(valid_samples)} samples: {', '.join(valid_samples)}")
    
    # Check JEC variations
    print("\nJEC Variation Status:")
    for sample_name in valid_samples:
        events = all_events[sample_name]
        print(f"  {sample_name}: JEC_up={events.has_jec_up}, JEC_down={events.has_jec_down}")
    
    # Make plots
    print("\n" + "="*80)
    print("Generating Plots")
    print("="*80)
    
    plot_pt_distributions(all_events)
    plot_response(all_events)
    
    print("\n" + "="*80)
    print("Done!")
    print("="*80)


if __name__ == '__main__':
    main()
