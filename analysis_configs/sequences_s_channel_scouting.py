import numpy as np
import awkward as ak

from skimmer import skimmer_utils
from utils.variables_computation import event_variables as event_vars
from utils.variables_computation import jet_variables as jet_vars
from utils.data.triggers import primary_dataset_triggers
from utils.tree_maker.triggers import trigger_table
from utils.awkward_array_utilities import as_type
import analysis_configs.physics_objects_nano_aod as physicsObj
from analysis_configs import objects_definition_s_channel_scouting as obj

from utils.Logger import *


def __unflatten_to_event_level(jet_variable,n_jets):
    arrays = [ak.unflatten(jet_variable, n_jets, axis=0)]

    if len(arrays) == 1:
        return arrays[0]
    else:
        return arrays

def apply_good_ak8_jet_filter(events):
    events = add_branches_for_ak8_jet_id(events)
    events = add_jet_id_branch(events)
    ak8_jets = ak.zip(
        {
            "pt": events.FatJet_pt,
            "mass": events.FatJet_mass,
            "eta": events.FatJet_eta,
            "phi": events.FatJet_phi,
            "id": events.FatJet_jetId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    
    analysis_jets = ak8_jets[obj.is_analysis_ak8_jet(ak8_jets)]
    good_ak8_jets_filter = ak.all(obj.is_good_ak8_jet(analysis_jets), axis=1)
    events = events[good_ak8_jets_filter]
    return events


def add_branches_for_ak8_jet_id(events):
        
        jet_collection_name = "FatJet"
        constituent_collection_name = "PFCands"
        jets = physicsObj.get_jets(events, jet_collection_name)
        jet_pf_cands = physicsObj.get_jet_pf_cands(events, jet_collection_name, constituent_collection_name)
        jet_pf_cands_per_jet, _ = jet_vars.make_constituents_per_jet(jet_pf_cands, n_jets=ak.count(events[jet_collection_name+"_pt"], axis=-1))
        flat_jets = ak.flatten(jets)

        variables_fractions = ["chHEF", "neHEF", "electron_energy_fraction", "muon_energy_fraction", "photon_energy_fraction"]

        for variable_name in variables_fractions:
            function = getattr(jet_vars, f"calculate_{variable_name}")
            variable = function(jet_pf_cands_per_jet, flat_jets)
            branch_name = f"FatJet_{variable_name.replace('_','')}"
            events[branch_name] = __unflatten_to_event_level(variable, n_jets=ak.count(events[jet_collection_name+"_pt"], axis=-1))

        variables_multiplicities = ["electron_multiplicity", "muon_multiplicity", "chargedhadron_multiplicity"]

        for variable_name in variables_multiplicities:
            function = getattr(jet_vars, f"calculate_{variable_name}")
            variable = function(jet_pf_cands_per_jet)
            branch_name = f"{jet_collection_name}_{variable_name.replace('_','')}"
            events[branch_name] = __unflatten_to_event_level(variable, n_jets=ak.count(events[jet_collection_name+"_pt"], axis=-1))

        #total PFCands multiplicity per jet
        function = getattr(jet_vars, "calculate_multiplicity")
        variable = function(jet_pf_cands_per_jet)
        branch_name = f"{jet_collection_name}_multiplicity"
        events[branch_name] = __unflatten_to_event_level(variable, n_jets=ak.count(events[jet_collection_name+"_pt"], axis=-1))

        return events

def add_jet_id_branch(events):
    #CZZ: here custom implementation of id function from offline analysis, to be substituted with variable in ntuplizer
    # passID
    CHM = events.FatJet_chargedhadronmultiplicity +events.FatJet_electronmultiplicity + events.FatJet_muonmultiplicity

    passID = ((abs(events.FatJet_eta)<=2.6) 
                   & (CHM>0) 
                   & (events.FatJet_chHEF>0) 
                   & (events.FatJet_neHEF<0.9) 
                   & (events.FatJet_multiplicity>1)
                   & (events.FatJet_photonenergyfraction<0.9)
                   & (events.FatJet_muonenergyfraction<0.8) 
                   & (events.FatJet_electronenergyfraction<0.8))  

    #convert passID to 0 and 1
    passID = as_type(passID,int)

    events["FatJet_jetId"] = passID

    return events

def add_good_ak8_jet_branch(events):
    ak8_jets = ak.zip(
        {
            "pt": events.FatJet_pt,
            "mass": events.FatJet_mass,
            "eta": events.FatJet_eta,
            "phi": events.FatJet_phi,
            "id": events.FatJet_jetId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    
    is_good_analysis_ak8_jet = (
        obj.is_analysis_ak8_jet(ak8_jets)
        & obj.is_good_ak8_jet(ak8_jets)
    )

    #add new branch to the events
    events["FatJet_isGood"] = is_good_analysis_ak8_jet

    return events


def _build_scouting_lepton_collections(events):
    required_fields = [
        "Electron_pt",
        "Electron_eta",
        "Electron_combinedMiniIso",
        "Muon_pt",
        "Muon_eta",
        "Muon_combinedMiniIso",
    ]

    if any(field not in events.fields for field in required_fields):
        return None, None

    electrons = ak.zip(
        {
            "pt": events["Electron_pt"],
            "eta": events["Electron_eta"],
            "iso": events["Electron_combinedMiniIso"],
        }
    )

    muons = ak.zip(
        {
            "pt": events["Muon_pt"],
            "eta": events["Muon_eta"],
            "iso": events["Muon_combinedMiniIso"],
        }
    )

    return electrons, muons

def __get_number_of_veto_leptons(
        events,
        electron_extra_condition=None,
        muon_extra_condition=None,
    ):

    electrons, muons = _build_scouting_lepton_collections(events)
    if electrons is None or muons is None:
        return ak.zeros_like(ak.num(events.FatJet_pt, axis=1), dtype=np.int64)

    electron_condition = obj.is_veto_electron(electrons)
    if electron_extra_condition is not None:
        electron_condition = electron_condition & electron_extra_condition

    muon_condition = obj.is_veto_muon(muons)
    if muon_extra_condition is not None:
        muon_condition = muon_condition & muon_extra_condition

    veto_electrons = electrons[electron_condition]
    veto_muons = muons[muon_condition]
    n_veto_electrons = ak.count(veto_electrons.pt, axis=1)
    n_veto_muons = ak.count(veto_muons.pt, axis=1)
    n_veto_leptons = n_veto_electrons + n_veto_muons

    return n_veto_leptons


def add_n_lepton_veto_branch(events):
    n_veto_leptons = __get_number_of_veto_leptons(
        events,
        electron_extra_condition=None,
        muon_extra_condition=None,
    )

    events = ak.with_field(
        events,
        n_veto_leptons,
        "nVetoLeptons",
    )

    return events
 

def apply_isolated_lepton_veto(
        events,
        electron_extra_condition=None,
        muon_extra_condition=None,
        revert=False,
    ):
    n_veto_leptons = __get_number_of_veto_leptons(
        events,
        electron_extra_condition=electron_extra_condition,
        muon_extra_condition=muon_extra_condition,
    )

    if revert:
        filter = (n_veto_leptons != 0)
    else:
        filter = (n_veto_leptons == 0)
    events = events[filter]

    return events

def add_analysis_branches(events):

    # Event variables
    good_jets_ak8_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.FatJet_pt[events.FatJet_isGood],
        eta=events.FatJet_eta[events.FatJet_isGood],
        phi=events.FatJet_phi[events.FatJet_isGood],
        mass=events.FatJet_mass[events.FatJet_isGood],
    )
    met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.ScoutMET_pt,
            phi=events.ScoutMET_phi,
    )

    #add rt variable
    mt = event_vars.calculate_transverse_mass(good_jets_ak8_lv, met)
    rt = met.pt / mt
    events["RTFatJet"] = rt

    #add mt variable
    events["MT01FatJetMET"] = mt

    #add deltaeta the two leading jets
    events["DeltaEtaJ0J1FatJet"] = event_vars.calculate_delta_eta(good_jets_ak8_lv)
    
    #add minimum delta phi between the MET and the two leading jets
    events["DeltaPhiMinFatJetMET"] = event_vars.calculate_delta_phi_min(good_jets_ak8_lv, met)
    
    return events


def has_dark_quark_info(events):
    return "MatrixElementGenParticle_pt" in events.fields


def add_dark_quark_matching(events):
    """
    performs delta R matching between AK8 jets and dark quarks and saves the indices of the matched jets in a new branch
    """

    jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.FatJet_pt,
        eta=events.FatJet_eta,
        phi=events.FatJet_phi,
        mass=events.FatJet_mass,
    )

    dark_quarks =  skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.MatrixElementGenParticle_pt,
        eta=events.MatrixElementGenParticle_eta,
        phi=events.MatrixElementGenParticle_phi,
        mass=events.MatrixElementGenParticle_mass,
    )

    # delta R matching between jets and dark quarks
    dR1 = jets.delta_r(dark_quarks[:,0]) <= 0.8
    dR2 = jets.delta_r(dark_quarks[:,1]) <= 0.8
    matched_mask = dR1 | dR2

    #add new branch to the event tree
    events["FatJet_DarkQuarkMatched"] = ak.values_astype(matched_mask, "uint32")
    
    return events


def remove_collections(events):
    #remove branch called genModel, hltResultName
    list_branches_to_remove = ["genModel", "hltResultName"]
    events = events[[key for key in events.fields if key not in list_branches_to_remove]]
    return events
