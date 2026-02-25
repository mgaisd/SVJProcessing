#!/usr/bin/env python3
"""
Create Jet Energy Correction maps for scouting jets.

This script matches scouting jets (Jet, FatJet) to offline jets (OfflineJet, OfflineFatJet)
and derives residual corrections as a function of jet pT and eta.

The correction map shows the ratio: pT(offline) / pT(scouting)

Performance optimizations:
- Vectorized jet matching using awkward arrays (10-100x faster)
- Single-pass file processing for both Jet and FatJet (2x I/O speedup)
- Numpy broadcasting for efficient ΔR calculations

Author: Created for SVJ scouting analysis
Date: February 2026
"""

import argparse
import os
import pickle
import json
from pathlib import Path
import subprocess

import numpy as np
import awkward as ak
import uproot
import matplotlib.pyplot as plt
import mplhep as hep

# Set CMS style for plots
plt.style.use(hep.style.CMS)


def get_file_list(dataset_path, use_xrootd=True):
    """
    Get list of ROOT files in the dataset directory.
    
    Args:
        dataset_path: Path to dataset (with or without xrootd redirector)
        use_xrootd: Whether to use xrootd (True) or local filesystem (False)
    
    Returns:
        List of file paths
    """
    if use_xrootd and not dataset_path.startswith("root://"):
        # Add default redirector
        redirector = "root://cmsdcache-kit-disk.gridka.de:1094/"
        full_path = f"{redirector}{dataset_path}"
    else:
        full_path = dataset_path
    
    if use_xrootd:
        # Extract server and path
        if full_path.startswith("root://"):
            parts = full_path.replace("root://", "").split("/", 1)
            server = parts[0]
            path = "/" + parts[1] if len(parts) > 1 else "/"
            
            # Use xrdfs to list files
            cmd = f"xrdfs {server} ls {path}"
            try:
                result = subprocess.run(cmd.split(), capture_output=True, text=True, timeout=30)
                files = [f"root://{server}/{line.strip()}" 
                         for line in result.stdout.split('\n') 
                         if line.strip().endswith('.root')]
                return files
            except Exception as e:
                print(f"Warning: xrdfs failed, trying glob: {e}")
                return []
    else:
        # Use local filesystem
        import glob
        files = sorted(glob.glob(os.path.join(dataset_path, "*.root")))
        return files


def match_jets_by_dr_vectorized(scout_pt, scout_eta, scout_phi, 
                                 offline_pt, offline_eta, offline_phi, 
                                 dr_threshold=0.1, max_jets=None, chunk_size=10000, verbose=True):
    """
    Match scouting jets to offline jets using fully vectorized awkward array operations.
    No Python loops - all operations are compiled C++ code for maximum speed.
    
    For each scouting jet, finds the closest offline jet within dr_threshold.
    
    Args:
        scout_pt, scout_eta, scout_phi: Scouting jet kinematics (awkward arrays)
        offline_pt, offline_eta, offline_phi: Offline jet kinematics (awkward arrays)
        dr_threshold: Maximum deltaR for matching
        max_jets: Maximum number of jets per event to consider (None = all jets)
        chunk_size: Process events in chunks to manage memory (default: 10000)
        verbose: Print matching efficiency (default: True)
    
    Returns:
        Tuple of (scouting_pt, scouting_eta, offline_pt_matched, n_matched, n_total)
        First three as flattened 1D numpy arrays with matched jets only
        Last two are statistics: number matched and total scouts
    """
    # DEBUG: Check input pT distributions
    scout_pt_flat_debug = ak.flatten(scout_pt)
    offline_pt_flat_debug = ak.flatten(offline_pt)
    if len(scout_pt_flat_debug) > 0:
        print(f"    DEBUG - Matching function input scout_pt: min={ak.min(scout_pt_flat_debug):.2f}, max={ak.max(scout_pt_flat_debug):.2f}, count={len(scout_pt_flat_debug)}")
    if len(offline_pt_flat_debug) > 0:
        print(f"    DEBUG - Matching function input offline_pt: min={ak.min(offline_pt_flat_debug):.2f}, max={ak.max(offline_pt_flat_debug):.2f}, count={len(offline_pt_flat_debug)}")
    
    # Apply max_jets cut if specified
    if max_jets is not None:
        scout_pt = scout_pt[:, :max_jets]
        scout_eta = scout_eta[:, :max_jets]
        scout_phi = scout_phi[:, :max_jets]
        offline_pt = offline_pt[:, :max_jets]
        offline_eta = offline_eta[:, :max_jets]
        offline_phi = offline_phi[:, :max_jets]
    
    n_events = len(scout_pt)
    n_chunks = (n_events + chunk_size - 1) // chunk_size
    
    all_scout_pt = []
    all_scout_eta = []
    all_offline_pt = []
    
    n_scout_total = 0
    n_matched = 0
    
    # Process in chunks to manage memory
    for chunk_idx in range(n_chunks):
        start_idx = chunk_idx * chunk_size
        end_idx = min((chunk_idx + 1) * chunk_size, n_events)
        
        # Get chunk
        scout_pt_chunk = scout_pt[start_idx:end_idx]
        scout_eta_chunk = scout_eta[start_idx:end_idx]
        scout_phi_chunk = scout_phi[start_idx:end_idx]
        offline_pt_chunk = offline_pt[start_idx:end_idx]
        offline_eta_chunk = offline_eta[start_idx:end_idx]
        offline_phi_chunk = offline_phi[start_idx:end_idx]
        
        # Compute deltaR for all scout-offline pairs within each event
        # Using broadcasting: (n_events, n_scout, 1) vs (n_events, 1, n_offline)
        deta = scout_eta_chunk[:, :, np.newaxis] - offline_eta_chunk[:, np.newaxis, :]
        dphi = scout_phi_chunk[:, :, np.newaxis] - offline_phi_chunk[:, np.newaxis, :]
        
        # Normalize phi to [-pi, pi] using awkward operations
        dphi = (dphi + np.pi) % (2 * np.pi) - np.pi
        
        # Compute deltaR: sqrt(deta^2 + dphi^2)
        dr = np.sqrt(deta**2 + dphi**2)
        
        # For each scout jet, find the offline jet with minimum deltaR
        # argmin along axis=2 (offline jets dimension)
        best_match_idx = ak.argmin(dr, axis=2, keepdims=False)
        min_dr = ak.min(dr, axis=2)
        
        # Apply deltaR threshold - creates boolean mask
        good_match = min_dr < dr_threshold
        
        # Get matched offline jets using the best match indices
        # Convert to global indices by flattening
        offline_pt_flat = ak.flatten(offline_pt_chunk)
        counts_offline = ak.num(offline_pt_chunk, axis=1)
        
        # Compute cumulative offsets for global indexing (vectorized)
        offsets = np.cumsum([0] + list(counts_offline))[:-1]  # [0, counts[0], counts[0]+counts[1], ...]
        offsets_ak = ak.unflatten(offsets, 1)  # Make it (n_events, 1) for broadcasting
        
        # Convert relative indices to global indices
        best_match_global = best_match_idx + offsets_ak
        
        # Flatten and index
        best_match_global_flat = ak.flatten(best_match_global)
        matched_offline_pt_flat = offline_pt_flat[best_match_global_flat]
        
        # Unflatten back to match scout structure
        counts_scout = ak.num(best_match_idx, axis=1)
        matched_offline_pt = ak.unflatten(matched_offline_pt_flat, counts_scout)
        
        # Apply good match mask to get only valid matches
        scout_pt_matched = scout_pt_chunk[good_match]
        scout_eta_matched = scout_eta_chunk[good_match]
        offline_pt_matched = matched_offline_pt[good_match]
        
        # Convert to numpy and accumulate
        all_scout_pt.append(ak.to_numpy(ak.flatten(scout_pt_matched)))
        all_scout_eta.append(ak.to_numpy(ak.flatten(scout_eta_matched)))
        all_offline_pt.append(ak.to_numpy(ak.flatten(offline_pt_matched)))
        
        # Track statistics
        n_scout_total += ak.sum(ak.num(scout_pt_chunk, axis=1))
        n_matched += ak.sum(ak.num(scout_pt_matched, axis=1))
    
    if verbose and n_scout_total > 0:
        eff = 100.0 * n_matched / n_scout_total
        print(f"    Matched {n_matched}/{n_scout_total} jets ({eff:.1f}%)")
    
    # Concatenate all chunks
    final_scout_pt = np.concatenate(all_scout_pt) if all_scout_pt else np.array([])
    final_scout_eta = np.concatenate(all_scout_eta) if all_scout_eta else np.array([])
    final_offline_pt = np.concatenate(all_offline_pt) if all_offline_pt else np.array([])
    
    # DEBUG: Check output pT distributions
    if len(final_scout_pt) > 0:
        print(f"    DEBUG - Matching function output scout_pt: min={np.min(final_scout_pt):.2f}, max={np.max(final_scout_pt):.2f}, count={len(final_scout_pt)}")
    if len(final_offline_pt) > 0:
        print(f"    DEBUG - Matching function output offline_pt: min={np.min(final_offline_pt):.2f}, max={np.max(final_offline_pt):.2f}, count={len(final_offline_pt)}")
    
    return (final_scout_pt, final_scout_eta, final_offline_pt, n_matched, n_scout_total)


def create_histogram_stats(pt_bins, eta_bins):
    """
    Create empty histogram statistics for incremental accumulation.
    
    Args:
        pt_bins: Bin edges for pT
        eta_bins: Bin edges for eta
    
    Returns:
        Dictionary with arrays for tracking weighted statistics per bin
    """
    n_pt_bins = len(pt_bins) - 1
    n_eta_bins = len(eta_bins) - 1
    
    return {
        'sum_ratios': np.zeros((n_pt_bins, n_eta_bins)),
        'sum_ratios_squared': np.zeros((n_pt_bins, n_eta_bins)),
        'counts': np.zeros((n_pt_bins, n_eta_bins)),
        'pt_bins': pt_bins,
        'eta_bins': eta_bins,
        # Quality control tracking
        'n_total_jets': 0,
        'n_rejected_negative_pt': 0,
        'n_rejected_outliers': 0
    }


def update_histogram_stats(stats, scout_pt, scout_eta, offline_pt, jet_type="Unknown", filename=""):
    """
    Update histogram statistics with new matched jets (memory-efficient).
    Instead of storing all jets, accumulates sums per bin (unweighted).
    
    Args:
        stats: Dictionary from create_histogram_stats()
        scout_pt: Array of scouting jet pT for this batch
        scout_eta: Array of scouting jet eta
        offline_pt: Array of matched offline jet pT
        jet_type: Jet collection type (e.g., "Jet", "FatJet") for debug output
        filename: Name of the file being processed (for debug output)
    
    Returns:
        Tuple of (n_rejected_negative_pt, n_rejected_outliers) for this batch
    """
    if len(scout_pt) == 0:
        return (0, 0)
    
    pt_bins = stats['pt_bins']
    eta_bins = stats['eta_bins']
    n_pt_bins = len(pt_bins) - 1
    n_eta_bins = len(eta_bins) - 1
    
    # Data quality checks
    n_total = len(scout_pt)
    stats['n_total_jets'] += n_total
    
    valid_mask = (scout_pt > 0) & (offline_pt > 0) & np.isfinite(scout_pt) & np.isfinite(offline_pt)
    n_invalid_pt = n_total - np.sum(valid_mask)
    
    if n_invalid_pt > 0:
        stats['n_rejected_negative_pt'] += n_invalid_pt
        
        # DEBUG: Print details about invalid jets
        invalid_scout_pt = scout_pt[~valid_mask]
        invalid_offline_pt = offline_pt[~valid_mask]
        jet_label = "AK4" if jet_type == "Jet" else "AK8" if jet_type == "FatJet" else jet_type
        print(f"\n    DEBUG - Invalid jets detected in {jet_label} [{filename}]:")
        print(f"      Scout pT <= 0: {np.sum(scout_pt <= 0)}, NaN/Inf: {np.sum(~np.isfinite(scout_pt))}")
        print(f"      Offline pT <= 0: {np.sum(offline_pt <= 0)}, NaN/Inf: {np.sum(~np.isfinite(offline_pt))}")
        print(f"      Example rejected jets (first {min(3, len(invalid_scout_pt))})")
        for i in range(min(3, len(invalid_scout_pt))):
            print(f"        Jet {i+1}: scout_pt={invalid_scout_pt[i]:.6f} GeV, offline_pt={invalid_offline_pt[i]:.6e} GeV")
        
        scout_pt = scout_pt[valid_mask]
        scout_eta = scout_eta[valid_mask]
        offline_pt = offline_pt[valid_mask]
    
    if len(scout_pt) == 0:
        return (n_invalid_pt, 0)
    
    # Calculate ratio for all jets
    ratio = offline_pt / scout_pt
    
    # Outlier rejection: reject unphysical correction factors
    # Reasonable range: 0.1 < correction < 10.0
    outlier_mask = (ratio > 0.1) & (ratio < 10.0)
    n_outliers = len(ratio) - np.sum(outlier_mask)
    
    if n_outliers > 0:
        stats['n_rejected_outliers'] += n_outliers
        
        # DEBUG: Print details about outlier rejections
        outlier_ratios = ratio[~outlier_mask]
        outlier_scout_pt = scout_pt[~outlier_mask]
        outlier_offline_pt = offline_pt[~outlier_mask]
        
        jet_label = "AK4" if jet_type == "Jet" else "AK8" if jet_type == "FatJet" else jet_type
        print(f"\n    DEBUG - Outlier ratios detected in {jet_label} [{filename}]:")
        print(f"      Count: {n_outliers}")
        print(f"      Ratio range: [{np.min(outlier_ratios):.6f}, {np.max(outlier_ratios):.6f}]")
        print(f"      Example rejected jets (first {min(5, len(outlier_ratios))})")
        for i in range(min(5, len(outlier_ratios))):
            print(f"        Jet {i+1}: scout_pt={outlier_scout_pt[i]:.2f} GeV, offline_pt={outlier_offline_pt[i]:.2f} GeV, ratio={outlier_ratios[i]:.6f}")
        
        scout_pt = scout_pt[outlier_mask]
        scout_eta = scout_eta[outlier_mask]
        offline_pt = offline_pt[outlier_mask]
        ratio = ratio[outlier_mask]
    
    if len(scout_pt) == 0:
        return (n_invalid_pt, n_outliers)
    
    # Loop over bins and accumulate weighted statistics
    for i_pt in range(n_pt_bins):
        for i_eta in range(n_eta_bins):
            # Create mask for this bin
            mask = ((scout_pt >= pt_bins[i_pt]) & 
                   (scout_pt < pt_bins[i_pt + 1]) &
                   (scout_eta >= eta_bins[i_eta]) & 
                   (scout_eta < eta_bins[i_eta + 1]))
            
            n_jets = np.sum(mask)
            
            if n_jets > 0:
                ratios_in_bin = ratio[mask]
                # Accumulate sums for this bin
                stats['sum_ratios'][i_pt, i_eta] += np.sum(ratios_in_bin)
                stats['sum_ratios_squared'][i_pt, i_eta] += np.sum(ratios_in_bin**2)
                stats['counts'][i_pt, i_eta] += n_jets
    
    return (n_invalid_pt, n_outliers)


def finalize_correction_map(stats):
    """
    Calculate final correction map from accumulated histogram statistics.
    
    Args:
        stats: Dictionary from create_histogram_stats() after updates
    
    Returns:
        Dictionary with:
            'corrections': 2D array of correction factors
            'uncertainties': 2D array of weighted standard deviations
            'counts': 2D array of effective number of jets per bin
            'pt_bins': pT bin edges
            'eta_bins': eta bin edges
    """
    n_pt_bins = stats['sum_ratios'].shape[0]
    n_eta_bins = stats['sum_ratios'].shape[1]
    
    corrections = np.ones((n_pt_bins, n_eta_bins))
    uncertainties = np.zeros((n_pt_bins, n_eta_bins))
    stat_uncertainties = np.zeros((n_pt_bins, n_eta_bins))
    jet_counts = np.zeros((n_pt_bins, n_eta_bins))  # Will copy from stats['counts']
    
    # Calculate mean and std for each bin
    n_problematic_bins = 0
    for i_pt in range(n_pt_bins):
        for i_eta in range(n_eta_bins):
            n = stats['counts'][i_pt, i_eta]
            
            if n > 0:
                # Mean of ratios
                mean = stats['sum_ratios'][i_pt, i_eta] / n
                
                # Sanity check: corrections should be positive and reasonable
                if mean <= 0 or mean < 0.1 or mean > 10.0:
                    n_problematic_bins += 1
                    pt_low = stats['pt_bins'][i_pt]
                    pt_high = stats['pt_bins'][i_pt+1]
                    eta_low = stats['eta_bins'][i_eta]
                    eta_high = stats['eta_bins'][i_eta+1]
                    print(f"  WARNING: Problematic bin at pT=[{pt_low:.0f},{pt_high:.0f}], eta=[{eta_low:.2f},{eta_high:.2f}]: "
                          f"correction={mean:.3f}, n_jets={int(n)}")
                    # Set to 1.0 (no correction) for problematic bins
                    mean = 1.0
                
                corrections[i_pt, i_eta] = mean
                
                # Standard deviation of ratios: sqrt(E[X²] - E[X]²)
                mean_of_squares = stats['sum_ratios_squared'][i_pt, i_eta] / n
                variance = mean_of_squares - mean**2
                uncertainties[i_pt, i_eta] = np.sqrt(max(0, variance))  # Resolution/response spread (JEC systematic)
                
                # Statistical uncertainty on the mean correction: std/sqrt(N_jets)
                if n > 1:
                    stat_uncertainties[i_pt, i_eta] = uncertainties[i_pt, i_eta] / np.sqrt(n)
                
                # Store jet count
                jet_counts[i_pt, i_eta] = n
    
    if n_problematic_bins > 0:
        print(f"  WARNING: Set {n_problematic_bins} problematic bins to correction=1.0")
    
    # Print quality control summary
    n_total = stats['n_total_jets']
    n_neg = stats['n_rejected_negative_pt']
    n_outliers = stats['n_rejected_outliers']
    n_accepted = n_total - n_neg - n_outliers
    
    if n_total > 0:
        print(f"\nData Quality Summary:")
        print(f"  Total matched jets: {n_total}")
        print(f"  Rejected (non-positive pT): {n_neg} ({100.0*n_neg/n_total:.2f}%)")
        print(f"  Rejected (outlier ratios): {n_outliers} ({100.0*n_outliers/n_total:.2f}%)")
        print(f"  Accepted jets: {n_accepted} ({100.0*n_accepted/n_total:.2f}%)")
    
    return {
        'corrections': corrections,
        'uncertainties': uncertainties,
        'stat_uncertainties': stat_uncertainties,
        'counts': jet_counts,
        'pt_bins': stats['pt_bins'],
        'eta_bins': stats['eta_bins']
    }


def plot_correction_map(correction_map, jet_type, output_dir):
    """
    Create 2D visualization of binned correction map (CMS style).
    Shows only the average correction values per bin.
    
    Args:
        correction_map: Dictionary from calculate_correction_map()
        jet_type: String identifier for jet type (e.g., "AK4", "AK8")
        output_dir: Directory to save plots
    """
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract data
    pt_bins = correction_map['pt_bins']
    eta_bins = correction_map['eta_bins']
    corrections = correction_map['corrections']
    counts = correction_map['counts']
    
    # Mask bins with insufficient statistics
    corrections_masked = np.ma.masked_where(counts < 10, corrections)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot 2D binned correction map using pcolormesh
    mesh = ax.pcolormesh(pt_bins, eta_bins, corrections_masked.T,
                         cmap='viridis', vmin=0.5, vmax=1.5,
                         shading='flat', edgecolors='none')
    
    # Labels and styling
    ax.set_xlabel("Jet $p_T$ [GeV]", fontsize=18)
    ax.set_ylabel("Jet $\\eta$", fontsize=18)
    ax.set_xlim(pt_bins[0], pt_bins[-1])
    ax.set_ylim(eta_bins[0], eta_bins[-1])
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    # Colorbar
    cbar = plt.colorbar(mesh, ax=ax, label="Correction Factor (Offline/Scouting)")
    cbar.set_label("Correction Factor", fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    
    # CMS label
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    # Set proper jet type label
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Save
    output_file = os.path.join(output_dir, f"scouting_correction_map_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot: {output_file}")
    plt.close()
    
    # Create statistics heatmap
    plot_statistics_map(correction_map, jet_type, output_dir)
    
    # Create statistical uncertainty heatmap
    plot_stat_uncertainty_map(correction_map, jet_type, output_dir)
    
    # Also create a profile plot: correction vs pt
    plot_correction_profiles(correction_map, jet_type, output_dir)


def plot_statistics_map(correction_map, jet_type, output_dir):
    """
    Create 2D visualization of statistics (number of jets) per bin.
    
    Args:
        correction_map: Dictionary from calculate_correction_map()
        jet_type: String identifier for jet type (e.g., "AK4", "AK8")
        output_dir: Directory to save plots
    """
    # Extract data
    pt_bins = correction_map['pt_bins']
    eta_bins = correction_map['eta_bins']
    counts = correction_map['counts']
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Use log scale for better visualization of large dynamic range
    counts_plot = np.ma.masked_where(counts == 0, counts)
    
    # Plot 2D statistics map using pcolormesh with log scale
    mesh = ax.pcolormesh(pt_bins, eta_bins, counts_plot.T,
                         cmap='YlOrRd', norm=plt.matplotlib.colors.LogNorm(vmin=1, vmax=counts.max()),
                         shading='flat', edgecolors='none')
    
    # Labels and styling
    ax.set_xlabel("Jet $p_T$ [GeV]", fontsize=18)
    ax.set_ylabel("Jet $\\eta$", fontsize=18)
    ax.set_xlim(pt_bins[0], pt_bins[-1])
    ax.set_ylim(eta_bins[0], eta_bins[-1])
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    # Colorbar
    cbar = plt.colorbar(mesh, ax=ax, label="Number of Jets per Bin")
    cbar.set_label("Number of Jets per Bin", fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    
    # CMS label
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    # Set proper jet type label
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Add text annotation for bins with >10 jets
    n_good_bins = np.sum(counts > 10)
    ax.text(0.05, 0.05, f"Bins with >10 jets: {n_good_bins}/{counts.size}", 
            transform=ax.transAxes,
            fontsize=14, verticalalignment='bottom', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Save
    output_file = os.path.join(output_dir, f"scouting_statistics_map_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot: {output_file}")
    plt.close()


def plot_stat_uncertainty_map(correction_map, jet_type, output_dir):
    """
    Create 2D visualization of propagated statistical uncertainties per bin.
    
    Args:
        correction_map: Dictionary from finalize_correction_map()
        jet_type: String identifier for jet type (e.g., "AK4", "AK8")
        output_dir: Directory to save plots
    """
    # Extract data
    pt_bins = correction_map['pt_bins']
    eta_bins = correction_map['eta_bins']
    stat_uncertainties = correction_map['stat_uncertainties']
    counts = correction_map['counts']
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Mask bins with insufficient statistics
    stat_unc_masked = np.ma.masked_where(counts < 10, stat_uncertainties)
    
    # Plot 2D statistical uncertainty map using pcolormesh
    mesh = ax.pcolormesh(pt_bins, eta_bins, stat_unc_masked.T,
                         cmap='plasma', vmin=0, vmax=np.percentile(stat_uncertainties[counts >= 10], 95) if np.any(counts >= 10) else 0.1,
                         shading='flat', edgecolors='none')
    
    # Labels and styling
    ax.set_xlabel("Jet $p_T$ [GeV]", fontsize=18)
    ax.set_ylabel("Jet $\\eta$", fontsize=18)
    ax.set_xlim(pt_bins[0], pt_bins[-1])
    ax.set_ylim(eta_bins[0], eta_bins[-1])
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    # Colorbar
    cbar = plt.colorbar(mesh, ax=ax, label="Statistical Uncertainty")
    cbar.set_label("Statistical Uncertainty on Correction", fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    
    # CMS label
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    # Set proper jet type label
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Add text annotation for bins with >10 jets
    n_good_bins = np.sum(counts > 10)
    ax.text(0.05, 0.05, f"Bins with >10 jets: {n_good_bins}/{counts.size}", 
            transform=ax.transAxes,
            fontsize=14, verticalalignment='bottom', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Save
    output_file = os.path.join(output_dir, f"scouting_stat_uncertainty_map_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot: {output_file}")
    plt.close()


def plot_correction_profiles(correction_map, jet_type, output_dir):
    """
    Create 1D profile plots showing corrections vs pT and vs eta.
    """
    pt_bins = correction_map['pt_bins']
    eta_bins = correction_map['eta_bins']
    corrections = correction_map['corrections']
    uncertainties = correction_map['uncertainties']
    counts = correction_map['counts']
    
    pt_centers = 0.5 * (pt_bins[:-1] + pt_bins[1:])
    eta_centers = 0.5 * (eta_bins[:-1] + eta_bins[1:])
    
    # Profile vs pT (averaged over eta)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for i_eta in range(len(eta_centers)):
        eta_label = f"{eta_bins[i_eta]:.1f} < η < {eta_bins[i_eta+1]:.1f}"
        
        # Extract corrections for this eta slice
        corr_vs_pt = corrections[:, i_eta]
        unc_vs_pt = uncertainties[:, i_eta]
        
        # Only plot bins with data
        mask = counts[:, i_eta] > 10
        if np.sum(mask) > 0:
            ax.errorbar(pt_centers[mask], corr_vs_pt[mask], 
                       yerr=unc_vs_pt[mask],
                       marker='o', linestyle='-', label=eta_label, alpha=0.7)
    
    ax.set_xlabel("Scouting Jet $p_T$ [GeV]", fontsize=18)
    ax.set_ylabel("Correction Factor (Offline/Scouting)", fontsize=18)
    ax.axhline(1.0, color='black', linestyle='--', alpha=0.5)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.legend(fontsize=12, ncol=2)
    ax.grid(alpha=0.3)
    
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    # Set proper jet type label
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    output_file = os.path.join(output_dir, f"scouting_correction_vs_pt_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot: {output_file}")
    plt.close()
    
    # Profile vs eta (averaged over pt)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Average over pT bins
    corr_vs_eta = np.nanmean(corrections, axis=0)
    
    # Calculate uncertainties, protecting against division by zero
    n_bins_with_data = np.sum(counts > 0, axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        unc_vs_eta = np.sqrt(np.nansum(uncertainties**2, axis=0)) / n_bins_with_data
        unc_vs_eta = np.where(n_bins_with_data > 0, unc_vs_eta, 0)  # Set invalid bins to 0
    
    mask = np.sum(counts, axis=0) > 10
    ax.errorbar(eta_centers[mask], corr_vs_eta[mask], 
               yerr=unc_vs_eta[mask],
               marker='o', linestyle='-', color='tab:blue')
    
    ax.set_xlabel("Scouting Jet $\\eta$", fontsize=18)
    ax.set_ylabel("Correction Factor (Offline/Scouting)", fontsize=18)
    ax.axhline(1.0, color='black', linestyle='--', alpha=0.5)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(alpha=0.3)
    
    hep.cms.label(loc=0, data=False, lumi=None, year=2017, ax=ax)
    
    # Set proper jet type label
    jet_label = "AK4 Jets" if jet_type == "Jet" else "AK8 Jets"
    ax.text(0.95, 0.95, jet_label, transform=ax.transAxes,
            fontsize=16, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    output_file = os.path.join(output_dir, f"scouting_correction_vs_eta_{jet_type}.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot: {output_file}")
    plt.close()


def count_events_in_dataset(files, use_xrootd=True):
    """
    Count total number of events in dataset files.
    
    Args:
        files: List of file paths
        use_xrootd: Whether files are on xrootd
    
    Returns:
        Total number of events
    """
    total_events = 0
    for file_path in files:
        try:
            with uproot.open(file_path) as f:
                tree = f['mmtree/Events']
                total_events += tree.num_entries
        except Exception as e:
            print(f"    Warning: Could not count events in {Path(file_path).name}: {e}")
    return total_events


def process_dataset(dataset_path, histogram_stats, jet_types=["Jet", "FatJet"], 
                    dr_threshold=0.1, max_files=None, use_xrootd=True, max_jets=None):
    """
    Process a QCD dataset and update histogram statistics incrementally.
    Memory-efficient: does not store all matched jets, only accumulates binned statistics.
    
    Args:
        dataset_path: Path to dataset directory
        histogram_stats: Dictionary of histogram stats per jet type (modified in-place)
        jet_types: List of jet types to process ["Jet", "FatJet"]
        dr_threshold: Maximum deltaR for matching
        max_files: Maximum number of files to process (None = all)
        use_xrootd: Whether to use xrootd
        max_jets: Maximum number of jets per event to consider (None = all)
    
    Returns:
        Number of events processed
    """
    print(f"\nProcessing dataset: {dataset_path}")
    print(f"  Jet types: {jet_types}")
    
    # Get file list
    files = get_file_list(dataset_path, use_xrootd=use_xrootd)
    if not files:
        print("  No files found!")
        return 0
    
    print(f"  Found {len(files)} files")
    
    if max_files is not None:
        files = files[:max_files]
        print(f"  Processing first {max_files} files")
    
    # Count total events for statistics
    n_events_processed = 0
    
    # Track matching statistics per jet type
    matching_stats = {jt: {'n_matched': 0, 'n_total': 0} for jt in jet_types}
    
    # Track if we've warned about missing branches
    warned_missing_branches = set()
    
    # Track files with rejections
    files_with_rejections = []
    
    # Process each file with simple progress indicator
    for i_file, file_path in enumerate(files, 1):
        # Update progress after every file
        pct = 100.0 * i_file / len(files)
        print(f"\r  Processing: {i_file}/{len(files)} ({pct:.1f}%)", end='', flush=True)
        
        try:
            file_has_rejections = False
            
            with uproot.open(file_path) as f:
                tree = f['mmtree/Events']
                branches = tree.keys()
                n_events_processed += tree.num_entries
                
                # Process each jet type
                for jet_type in jet_types:
                    # Define branch names
                    if jet_type == "Jet":
                        scout_prefix = "Jet"
                        offline_prefix = "OfflineJet"
                    elif jet_type == "FatJet":
                        scout_prefix = "FatJet"
                        offline_prefix = "OfflineFatJet"
                    else:
                        continue
                    
                    # Check required branches
                    required = [f"{scout_prefix}_pt", f"{scout_prefix}_eta", f"{scout_prefix}_phi",
                               f"{offline_prefix}_pt", f"{offline_prefix}_eta", f"{offline_prefix}_phi"]
                    
                    missing = [br for br in required if br not in branches]
                    if missing:
                        # Only warn once per jet type
                        if jet_type not in warned_missing_branches:
                            print(f"\n    Warning ({jet_type}): Missing branches: {missing}")
                            warned_missing_branches.add(jet_type)
                        continue
                    
                    # Load data
                    scout_pt = tree[f"{scout_prefix}_pt"].array(library="ak")
                    scout_eta = tree[f"{scout_prefix}_eta"].array(library="ak")
                    scout_phi = tree[f"{scout_prefix}_phi"].array(library="ak")
                    
                    offline_pt = tree[f"{offline_prefix}_pt"].array(library="ak")
                    offline_eta = tree[f"{offline_prefix}_eta"].array(library="ak")
                    offline_phi = tree[f"{offline_prefix}_phi"].array(library="ak")
                    
                    # DEBUG: Check for low pT AK8 jets in raw data
                    if jet_type == "FatJet":
                        scout_pt_flat = ak.flatten(scout_pt)
                        offline_pt_flat = ak.flatten(offline_pt)
                        
                        low_scout = scout_pt_flat < 100
                        low_offline = offline_pt_flat < 100
                        
                        n_low_scout = ak.sum(low_scout)
                        n_low_offline = ak.sum(low_offline)
                        
                        if n_low_scout > 0 or n_low_offline > 0:
                            print(f"\n    *** RAW DATA CHECK [{Path(file_path).name}] ***")
                            if n_low_scout > 0:
                                print(f"    Found {n_low_scout} FatJet (scout) with pT < 100 GeV in raw data!")
                                low_pts = scout_pt_flat[low_scout]
                                print(f"      pT range: [{ak.min(low_pts):.2f}, {ak.max(low_pts):.2f}] GeV")
                                print(f"      First 10: {ak.to_numpy(low_pts[:10])}")
                            if n_low_offline > 0:
                                print(f"    Found {n_low_offline} OfflineFatJet with pT < 100 GeV in raw data!")
                                low_pts_off = offline_pt_flat[low_offline]
                                print(f"      pT range: [{ak.min(low_pts_off):.2f}, {ak.max(low_pts_off):.2f}] GeV")
                                print(f"      First 10: {ak.to_numpy(low_pts_off[:10])}")
                        else:
                            print(f"\n    RAW DATA CHECK [{Path(file_path).name}]: No AK8 jets with pT < 100 GeV found (all jets have pT >= 100 GeV)")                    
                    # Apply pT cuts: AK4 jets must have pT > 15 GeV, AK8 jets > 150 GeV
                    if jet_type == "Jet":
                        # Cut on scout jets
                        scout_mask = scout_pt > 15
                        scout_pt = scout_pt[scout_mask]
                        scout_eta = scout_eta[scout_mask]
                        scout_phi = scout_phi[scout_mask]
                        
                        # Cut on offline jets
                        offline_mask = offline_pt > 15
                        offline_pt = offline_pt[offline_mask]
                        offline_eta = offline_eta[offline_mask]
                        offline_phi = offline_phi[offline_mask]
                    elif jet_type == "FatJet":
                        # Cut on scout jets
                        scout_mask = scout_pt > 150
                        scout_pt = scout_pt[scout_mask]
                        scout_eta = scout_eta[scout_mask]
                        scout_phi = scout_phi[scout_mask]
                        
                        # Cut on offline jets
                        offline_mask = offline_pt > 150
                        offline_pt = offline_pt[offline_mask]
                        offline_eta = offline_eta[offline_mask]
                        offline_phi = offline_phi[offline_mask]
                    
                    # Match jets using vectorized function
                    # For FatJet, apply max_jets limit; for Jet, use all jets (Type-1 MET)
                    jets_limit = max_jets if jet_type == "FatJet" else None
                    scout_pt_matched, scout_eta_matched, offline_pt_matched, n_matched, n_total = match_jets_by_dr_vectorized(
                        scout_pt, scout_eta, scout_phi,
                        offline_pt, offline_eta, offline_phi,
                        dr_threshold=dr_threshold,
                        max_jets=jets_limit,
                        verbose=False
                    )
                    
                    # Accumulate matching statistics
                    matching_stats[jet_type]['n_matched'] += n_matched
                    matching_stats[jet_type]['n_total'] += n_total
                    
                    # Update histogram statistics incrementally (memory-efficient!)
                    n_rejected_negative, n_rejected_outliers = update_histogram_stats(
                        histogram_stats[jet_type],
                        scout_pt_matched, scout_eta_matched, offline_pt_matched,
                        jet_type=jet_type,
                        filename=Path(file_path).name
                    )
                    
                    # Track if this file has rejections
                    if n_rejected_negative > 0 or n_rejected_outliers > 0:
                        file_has_rejections = True
            
            # Record file if it had rejections
            if file_has_rejections:
                files_with_rejections.append(Path(file_path).name)
                
        except Exception as e:
            print(f"\n    Error processing file {Path(file_path).name}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Clear progress line
    print()
    
    # Print summary
    print(f"  Processed {n_events_processed} events from {len(files)} files")
    
    # Print files with rejections
    if files_with_rejections:
        print(f"\n  Files with rejected jets ({len(files_with_rejections)} total):")
        for filename in files_with_rejections[:10]:  # Show first 10
            print(f"    - {filename}")
        if len(files_with_rejections) > 10:
            print(f"    ... and {len(files_with_rejections) - 10} more files")
    
    # Print matching statistics for this dataset
    for jet_type in jet_types:
        stats = matching_stats[jet_type]
        if stats['n_total'] > 0:
            eff = 100.0 * stats['n_matched'] / stats['n_total']
            print(f"  {jet_type}: Matched {stats['n_matched']}/{stats['n_total']} jets ({eff:.1f}%)")
    
    return n_events_processed


def save_corrections(correction_map, jet_type, year, output_path):
    """
    Save correction map in multiple formats.
    
    Args:
        correction_map: Dictionary from calculate_correction_map()
        jet_type: String identifier for jet type
        year: Data-taking year
        output_path: Base path for output files
    """
    # Save as pickle (coffea format)
    coffea_data = {
        f'ScoutingJetCorrection_{jet_type}_{year}': correction_map
    }
    
    coffea_file = output_path.replace('.coffea', f'_{jet_type}.coffea')
    with open(coffea_file, 'wb') as f:
        pickle.dump(coffea_data, f)
    print(f"Saved coffea file: {coffea_file}")
    
    # Save as JSON (human-readable)
    json_data = {
        'jet_type': jet_type,
        'year': year,
        'pt_bins': correction_map['pt_bins'].tolist(),
        'eta_bins': correction_map['eta_bins'].tolist(),
        'corrections': correction_map['corrections'].tolist(),
        'uncertainties': correction_map['uncertainties'].tolist(),
        'stat_uncertainties': correction_map['stat_uncertainties'].tolist(),
        'counts': correction_map['counts'].tolist()
    }
    
    json_file = output_path.replace('.coffea', f'_{jet_type}.json')
    with open(json_file, 'w') as f:
        json.dump(json_data, f, indent=2)
    print(f"Saved JSON file: {json_file}")


def main():
    """Main function to orchestrate correction map creation."""
    parser = argparse.ArgumentParser(
        description='Create scouting jet energy correction maps',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--year', type=str, default='2017',
                       help='Data-taking year')
    parser.add_argument('--jet-type', type=str, default='both',
                       choices=['Jet', 'FatJet', 'both'],
                       help='Which jet collection to process')
    parser.add_argument('--dr-threshold', type=float, default=0.1,
                       help='Maximum deltaR for jet matching')
    parser.add_argument('--n-files', type=int, default=None,
                       help='Maximum number of files per dataset (for testing)')
    # Real output directories (production use):
    # parser.add_argument('--output-dir', type=str, default='data',
    #                    help='Output directory for correction files')
    # parser.add_argument('--plot-dir', type=str, default='jec_plots',
    #                    help='Output directory for diagnostic plots')
    
    # Test output directory:
    parser.add_argument('--output-dir', type=str, default='test',
                       help='Output directory for correction files')
    parser.add_argument('--plot-dir', type=str, default='test',
                       help='Output directory for diagnostic plots')
    parser.add_argument('--no-plot', action='store_true',
                       help='Skip creating diagnostic plots')
    parser.add_argument('--local', action='store_true',
                       help='Use local filesystem instead of xrootd')
    parser.add_argument('--pt-bins', type=int, default=30,
                       help='Number of pT bins')
    parser.add_argument('--eta-bins', type=int, default=24,
                       help='Number of eta bins')
    parser.add_argument('--max-jets', type=int, default=None,
                       help='Maximum jets per event for FatJet only (e.g., 2 for leading two). Jet (AK4) always uses all jets for Type-1 MET.')
    
    args = parser.parse_args()
    
    print("="*80)
    print("Scouting Jet Energy Correction Map Creator")
    print("="*80)
    print(f"Year: {args.year}")
    print(f"Jet type: {args.jet_type}")
    print(f"deltaR threshold: {args.dr_threshold}")
    print(f"pT bins: {args.pt_bins} (Jet: 15-2000 GeV, FatJet: 150-2000 GeV)")
    print(f"eta range: Jet: unlimited, FatJet: -2.4 to 2.4 ({args.eta_bins} bins)")
    if args.max_jets is not None:
        print(f"Using leading {args.max_jets} FatJet(s) per event (Jet uses all jets for Type-1 MET)")
    else:
        print("Using all jets per event")
    
    # Define QCD datasets
    qcd_datasets = [
        f"/store/user/mgaisdor/SVJScouting_ntuples/MC_offline/{args.year}/QCD_HT300to500",
        f"/store/user/mgaisdor/SVJScouting_ntuples/MC_offline/{args.year}/QCD_HT500to700",
        f"/store/user/mgaisdor/SVJScouting_ntuples/MC_offline/{args.year}/QCD_HT700to1000",
        f"/store/user/mgaisdor/SVJScouting_ntuples/MC_offline/{args.year}/QCD_HT1000to1500",
        f"/store/user/mgaisdor/SVJScouting_ntuples/MC_offline/{args.year}/QCD_HT1500to2000",
        f"/store/user/mgaisdor/SVJScouting_ntuples/MC_offline/{args.year}/QCD_HT2000toInf",
    ]
    
    # Determine which jet types to process
    jet_types = ['Jet', 'FatJet'] if args.jet_type == 'both' else [args.jet_type]
    
    # Define pT and eta binning for each jet type upfront
    binning = {}
    for jet_type in jet_types:
        if jet_type == "Jet":
            pt_bins = np.linspace(15, 2000, args.pt_bins + 1)
            eta_bins = np.linspace(-5.0, 5.0, args.eta_bins + 1)
        else:  # FatJet
            pt_bins = np.linspace(150, 2000, args.pt_bins + 1)
            eta_bins = np.linspace(-2.4, 2.4, args.eta_bins + 1)
        binning[jet_type] = (pt_bins, eta_bins)
    
    # Initialize histogram statistics (memory-efficient approach)
    histogram_stats = {}
    for jet_type in jet_types:
        pt_bins, eta_bins = binning[jet_type]
        histogram_stats[jet_type] = create_histogram_stats(pt_bins, eta_bins)
    
    print("\n" + "="*80)
    print("Processing QCD datasets")
    print("="*80)
    
    # Process datasets and accumulate histogram statistics
    print("\nProcessing datasets...")
    for dataset_path in qcd_datasets:
        # Process dataset and update histogram statistics incrementally
        n_processed = process_dataset(
            dataset_path,
            histogram_stats,
            jet_types=jet_types,
            dr_threshold=args.dr_threshold,
            max_files=args.n_files,
            use_xrootd=not args.local,
            max_jets=args.max_jets
        )
        print(f"  Processed {n_processed} events from {Path(dataset_path).name}")
    
    # Now finalize correction maps for each jet type
    for jet_type in jet_types:
        print("\n" + "="*80)
        print(f"Finalizing correction map for {jet_type} collection")
        print("="*80)
        
        # Define pT binning info for printing
        pt_bins, eta_bins = binning[jet_type]
        if jet_type == "Jet":
            print(f"Using pT range: 15-2000 GeV ({args.pt_bins} bins)")
            print(f"Using eta range: -5.0 to 5.0 ({args.eta_bins} bins)")
        else:  # FatJet
            print(f"Using pT range: 150-2000 GeV ({args.pt_bins} bins)")
            print(f"Using eta range: -2.4 to 2.4 ({args.eta_bins} bins)")
        
        # Finalize correction map from accumulated statistics
        print("\nCalculating correction map...")
        correction_map = finalize_correction_map(histogram_stats[jet_type])
        
        # Calculate total jets
        total_jets = np.sum(correction_map['counts'])
        print(f"Total jets used: {total_jets:.0f}")
        
        # Print summary
        print("\nCorrection Map Summary:")
        print(f"  Total bins: {correction_map['corrections'].size}")
        print(f"  Bins with any data: {np.sum(correction_map['counts'] > 0)}")
        print(f"  Bins with sufficient statistics (>10 jets, shown in plots): {np.sum(correction_map['counts'] > 10)}")
        bins_with_stats = correction_map['counts'] > 10
        print(f"  Mean correction: {np.mean(correction_map['corrections'][bins_with_stats]):.3f}")
        print(f"  Correction range: [{np.min(correction_map['corrections'][bins_with_stats]):.3f}, "
              f"{np.max(correction_map['corrections'][bins_with_stats]):.3f}]")
        
        # Save corrections
        output_path = os.path.join(args.output_dir, f"scouting_jec_residuals_{args.year}.coffea")
        save_corrections(correction_map, jet_type, args.year, output_path)
        
        # Create plots
        if not args.no_plot:
            print("\nCreating diagnostic plots...")
            plot_correction_map(correction_map, jet_type, args.plot_dir)
    
    print("\n" + "="*80)
    print("Done!")
    print("="*80)


if __name__ == "__main__":
    main()
