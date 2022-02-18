from __future__ import annotations

from typing import Any

def lepton_channel_categorization(
    events: Any,
    iso_muons: str,
    iso_electrons: str,
    era: str,
    is_mc: bool,
    run_period: str | None = None,
    noniso_muons: str | None = None,
    noniso_electrons: str | None = None,
) -> Any:
    """Define event-level masks categorizing events into different lepton channels (ee_OS, emu_SS, mumu_OS,
    maybe 1 and 3-lepton channels, or (subsets) of the 1-lepton channels with an additional nonisolated lepton"""
    # might be better to separate channel categorization from triggering, but lots of overlap calculations:
    # define the different channels, e.g. ee (OS), emu (SS), etc.

    if not iso_muons.endswith("_"):
        iso_muons += "_"
    if not iso_electrons.endswith("_"):
        iso_electrons += "_"
    if noniso_muons and not noniso_muons.endswith("_"):
        noniso_muons += "_"
    if oniso_electrons and not noniso_electrons.endswith("_"):
        noniso_electrons += "_"
        
    # hold tag for our iso lepton multiplicities
    n_iso_e = f"n{iso_electrons[:-1]}"
    n_iso_mu = f"n{iso_muons[:-1]}"
    n_noniso_e = f"n{noniso_electrons[:-1]}" if noniso_electrons else "0"
    n_noniso_mu = f"n{noniso_muons[:-1]}" if noniso_muons else "0"
    
    # Get our lepton charge sums for OSDL/SSDL categorization
    events = events.Define("sum_charge_iso_e", f"Sum({iso_electrons}charge)")
    events = events.Define("sum_charge_iso_mu", f"Sum({iso_muons}charge)")
    events = events.Define("sum_charge_iso_lep", "sum_charge_iso_e + sum_charge_iso_mu")

    # Get our noniso lepton multiplicities and charges, if applicable
    if noniso_muons:
        events = events.Define("sum_charge_noniso_mu", f"Sum({noniso_muons}charge)")
    else:
        events = events.Define("sum_charge_noniso_mu", "return 0;")

    if noniso_electrons:
        events = events.Define("sum_charge_noniso_e", f"Sum({noniso_electrons}charge)")
    else:
        events = events.Define("sum_charge_noniso_e", "return 0;")
    events = events.Define("sum_charge_noniso_lep", "sum_charge_noniso_e + sum_charge_noniso_mu")

    events = events.Define("sum_charge_all", "sum_charge_iso_lep + sum_charge_noniso_lep")
    
    
    # Single lepton
    events = events.Define("iso_1e0mu", f"({n_iso_e} == 1) && ({n_iso_mu} == 0)")
    events = events.Define("iso_0e1mu", f"({n_iso_e} == 0) && ({n_iso_mu} == 1)")

    # Dilepton, OS + SS not differentiated yet
    events = events.Define("iso_2e0mu", f"({n_iso_e} == 2) && ({n_iso_mu} == 0)")
    events = events.Define("iso_1e1mu", f"({n_iso_e} == 1) && ({n_iso_mu} == 1)")
    events = events.Define("iso_0e2mu", f"({n_iso_e} == 0) && ({n_iso_mu} == 2)")

    # Trilepton
    events = events.Define("iso_3e0mu", f"({n_iso_e} == 3) && ({n_iso_mu} == 0)")
    events = events.Define("iso_2e1mu", f"({n_iso_e} == 2) && ({n_iso_mu} == 1)")
    events = events.Define("iso_1e2mu", f"({n_iso_e} == 1) && ({n_iso_mu} == 2)")
    events = events.Define("iso_0e3mu", f"({n_iso_e} == 0) && ({n_iso_mu} == 3)")

    # Quadralepton
    events = events.Define("iso_4e0mu", f"({n_iso_e} == 4) && ({n_iso_mu} == 0)")
    events = events.Define("iso_3e1mu", f"({n_iso_e} == 3) && ({n_iso_mu} == 1)")
    events = events.Define("iso_2e2mu", f"({n_iso_e} == 2) && ({n_iso_mu} == 2)")
    events = events.Define("iso_1e3mu", f"({n_iso_e} == 1) && ({n_iso_mu} == 3)")
    events = events.Define("iso_0e4mu", f"({n_iso_e} == 0) && ({n_iso_mu} == 4)")

    # 1 isolated, 1 non-isolated leptons... for QCD background estimations... orthogonal to each other but not above channels
    events = events.Define(
        "iso_1e0mu_noniso_1e0mu",
        f"({n_iso_e} == 1) && ({n_iso_mu} == 0) && ({n_noniso_e} == 1) && ({n_noniso_mu} == 0)"
    )
    events = events.Define(
        "iso_1e0mu_noniso_0e1mu",
        f"({n_iso_e} == 1) && ({n_iso_mu} == 0) && ({n_noniso_e} == 0) && ({n_noniso_mu} == 1)"
    )
    events = events.Define(
        "iso_0e1mu_noniso_1e0mu",
        f"({n_iso_e} == 0) && ({n_iso_mu} == 1) && ({n_noniso_e} == 1) && ({n_noniso_mu} == 0)"
    )
    events = events.Define(
        "iso_0e1mu_noniso_0e1mu",
        f"({n_iso_e} == 0) && ({n_iso_mu} == 1) && ({n_noniso_e} == 0) && ({n_noniso_mu} == 1)"
    )

    # Add sum of isolated and nonisolated charges
    events = events.Define("iso_sumc0", "(sum_charge_iso_lep == 0)")
    events = events.Define("noniso_sumc0", "(sum_charge_noniso_lep == 0)")
    events = events.Define("sumc0", "(sum_charge_all == 0)")

    # Add the main channels of interest... this set should be orthogonal with each other
    events = events.Define("channel_e", "(iso_1e0mu == true)")
    events = events.Define("channel_mu", "(iso_0e1mu == true)")
    events = events.Define("channel_ee_OS", "(iso_2e0mu == true && iso_sumc0 == true)")
    events = events.Define("channel_emu_OS", "(iso_1e1mu == true && iso_sumc0 == true)")
    events = events.Define("channel_mumu_OS", "(iso_0e2mu == true && iso_sumc0 == true)")
    events = events.Define("channel_ee_SS", "(iso_2e0mu == true && iso_sumc0 == false)")
    events = events.Define("channel_emu_SS", "(iso_1e1mu == true && iso_sumc0 == false)")
    events = events.Define("channel_mumu_SS", "(iso_0e2mu == true && iso_sumc0 == false)")

    # Add the inverted iso channels of interest, SUBSET OF 'e' and 'mu' channels!
    # isolated-lepton_nonisolated-lepton_dilepton-charge format
    events = events.Define("channel_e_nie_OS", "(iso_1e0mu_noniso_1e0mu == true, sumc0 == true)")
    events = events.Define("channel_e_nie_SS", "(iso_1e0mu_noniso_1e0mu == true, sumc0 == false)")
    events = events.Define("channel_e_nim_OS", "(iso_1e0mu_noniso_0e1mu == true, sumc0 == true)")
    events = events.Define("channel_e_nim_SS", "(iso_1e0mu_noniso_0e1mu == true, sumc0 == false)")
    events = events.Define("channel_mu_nie_OS", "(iso_0e1mu_noniso_1e0mu == true, sumc0 == true)")
    events = events.Define("channel_mu_nie_SS", "(iso_0e1mu_noniso_1e0mu == true, sumc0 == false)")
    events = events.Define("channel_mu_nim_OS", "(iso_0e1mu_noniso_0e1mu == true, sumc0 == true)")
    events = events.Define("channel_mu_nim_SS", "(iso_0e1mu_noniso_0e1mu == true, sumc0 == false)")

    return events


def dilepton_trigger_selection(
    events: Any,
    iso_muons: str,
    iso_electrons: str,
    era: str,
    is_mc: bool,
    run_period: str | None = None,
    data_stream: str | None = None,
    noniso_muons: str | None = None,
    noniso_electrons: str | None = None,
) -> Any:    
    """
    Takes as input an events dataframe, collection names for isolated leptons, year, run period, data stream, and 
    WIP note: Could be written to take advantage of DefinePerSample, but would be a significant departure from coffea"""
    if not iso_muons.endswith("_"):
        iso_muons += "_"
    if not iso_electrons.endswith("_"):
        iso_electrons += "_"
    if noniso_muons and not noniso_muons.endswith("_"):
        noniso_muons += "_"
    if noniso_electrons and not noniso_electrons.endswith("_"):
        noniso_electrons += "_"
        
    events = events.Define("first_iso_muon_pt", 
        f"{iso_muons}pt.at(0, 0.0)"
    )
    events = events.Define("second_iso_muon_pt", 
        f"{iso_muons}pt.at(1, 0.0)"
    )
    events = events.Define("first_iso_electron_pt", 
        f"{iso_electrons}pt.at(0, 0.0)"
    )
    events = events.Define("second_iso_electron_pt", 
        f"{iso_electrons}pt.at(1, 0.0)"
    )

    events = events.Define(
        "subtrigger_MET",
        """(HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned
            || HLT_PFMET200_HBHECleaned
            || HLT_PFMET200_NotCleaned
           ) && (MET_pt > 210)"""
    )
    if era == "2018":
        events = events.Define(
            "trig_emu_MuonEG1",
                """(HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == true)
                && (first_iso_electron_pt > 25)
                && (first_iso_muon_pt > 15)"""
        )
        events = events.Define(
            "trig_emu_MuonEG2",
                """(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == true)
                && (first_iso_electron_pt > 15)
                && (first_iso_muon_pt > 25)"""
        )
        events = events.Define(
            "trig_mumu_DoubleMuon",
                """(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 == true)
                && (first_iso_muon_pt > 25)
                && (second_iso_muon_pt > 15)"""
        )
        events = events.Define(
            "trig_ee_EGamma1",
                """(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL == true)
                && (first_iso_electron_pt > 25)
                && (second_iso_electron_pt > 15)"""
        )
        # backup channels not used in the current analysis, but gains back some events...
        events = events.Define(
            "trig_emu_SingleMuon",
                """(HLT_IsoMu24 == true)
                && (first_iso_electron_pt > 15)
                && (first_iso_muon_pt > 27)"""
        )
        events = events.Define(
            "trig_emu_EGamma",
                """(HLT_Ele32_WPTight_Gsf == true)
                && (first_iso_electron_pt > 35)
                && (first_iso_muon_pt > 15)"""
        )
        events = events.Define(
            "trig_mumu_SingleMuon",
                """(HLT_IsoMu24 == true)
                && (first_iso_muon_pt > 27)
                && (second_iso_muon_pt > 15)"""
        )
        events = events.Define(
            "trig_ee_EGamma2",
                """(HLT_Ele32_WPTight_Gsf == true)
                && (first_iso_electron_pt > 35)
                && (second_iso_electron_pt > 15)"""
        )
        events = events.Define(
            "trig_emu_MET",
                """subtrigger_MET
                && ((first_iso_electron_pt > 25) && (first_iso_muon_pt > 15))
                || ((first_iso_electron_pt > 15) && (first_iso_muon_pt > 25))"""
        )
        events = events.Define(
            "trig_mumu_MET",
            """(subtrigger_MET && (first_iso_muon_pt > 25) && (second_iso_muon_pt > 15))"""
        )
        events = events.Define(
            "trig_ee_MET",
                """subtrigger_MET
                && (first_iso_electron_pt > 25)
                && (second_iso_electron_pt > 15)"""
        )
        if is_mc:
            events = events.Define(
                "trig_emu", "(trig_emu_MuonEG1 || trig_emu_MuonEG2 || trig_emu_SingleMuon || trig_emu_EGamma)"
            )
            events = events.Define(
                "trig_mumu", "(trig_mumu_DoubleMuon || trig_mumu_SingleMuon"
            )
            events = events.Define(
                "trig_ee", "(trig_ee_EGamma1 || trig_ee_EGamma2)"
            )
        else:
            if data_stream == "MuonEG":
                # Select emu if they pass any triggers, we'll pick up the rest from the exclusive events in SingleMuon and EGamma
                events = events.Define(
                    "trig_emu", "(trig_emu_MuonEG1 || trig_emu_MuonEG2)"
                )# , "trig_emu_SingleMuon || trig_emu_EGamma"),
                events = events.Define("trig_mumu", "false")
                events = events.Define("trig_ee", "false")
            elif data_stream == "DoubleMuon":
                events = events.Define("trig_emu", "false")
                events = events.Define(
                    "trig_mumu", "trig_mumu_DoubleMuon"
                )# , "trig_mumu_SingleMuon"),
                events = events.Define("trig_ee", "false")
            elif data_stream == "SingleMuon":
                events = events.Define(
                    "trig_emu",
                    "(trig_emu_MuonEG1  == false && trig_emu_MuonEG2  == false && trig_emu_SingleMuon == true)"
                )
                events = events.Define(
                    "trig_mumu",
                    "trig_mumu_DoubleMuon  == false && trig_mumu_SingleMuon == true)"
                )
                events = events.Define("trig_ee", "false")
            elif data_stream == "EGamma":
                events = events.Define(
                    "trig_emu",
                    """(trig_emu_MuonEG1  == false && trig_emu_MuonEG2  == false &&
                        trig_emu_SingleMuon  == false && trig_emu_EGamma == true) """
                )
                events = events.Define("trig_mumu", "false")
                events = events.Define(
                    "trig_ee", "(trig_ee_EGamma1 || trig_ee_EGamma2)"
                )
            elif data_stream == "MET":
                raise NotImplementedError
            else:
                raise ValueError(
                    f"data_stream must be specified for is_mc=False (for event overlap removal), input: {data_stream}"
                )
            # Do e.g. trig_masks.require(trig_emu_MuonEG1  == false &&, trig_emu_MuonEG2  == false &&, trig_emu_SingleMuon=True)
            # for the SingleMuon datastream/backup trigger, so no double counting events in the MuonEG datastream+trigger
    elif era == "2017":
        events = events.Define(
            "trig_emu_MuonEG1",
            (
                """(HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == true)
                && (first_iso_electron_pt > 25)
                && (first_iso_muon_pt > 15)"""
        )
        events = events.Define(
            "trig_emu_MuonEG2",
            (
                """(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == true)
                && (first_iso_electron_pt > 15)
                && (first_iso_muon_pt > 25)"""
        )
        if not is_mc and run_period == "B":
            events = events.Define(
                "trig_mumu_DoubleMuon",
                (
                    """(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ == true)
                    && (first_iso_muon_pt > 25)
                    && (second_iso_muon_pt > 15)"""
            )
        else:
            events = events.Define(
                "trig_mumu_DoubleMuon",
                (
                    """(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 == true)
                    && (first_iso_muon_pt > 25)
                    && (second_iso_muon_pt > 15)"""
            )
        events = events.Define(
            "trig_ee_DoubleEG",
            (
                """(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL == true)
                && (first_iso_electron_pt > 25)
                && (second_iso_electron_pt > 15""")
        )

        events = events.Define(
            "trig_emu_SingleMuon",
            (
                """(HLT_IsoMu27 == true)
                && (first_iso_electron_pt > 15)
                && (first_iso_muon_pt > 30)"""
        )
        events = events.Define(
            "trig_emu_SingleElectron",
            (
                """(HLT_Ele35_WPTight_Gsf == true)
                && (first_iso_electron_pt > 38)
                && (first_iso_muon_pt > 15)"""
        )
        events = events.Define(
            "trig_mumu_SingleMuon",
            (
                """(HLT_IsoMu27 == true)
                && (first_iso_muon_pt > 30)
                && (second_iso_muon_pt > 15)"""
        )
        events = events.Define(
            "trig_ee_SingleElectron",
            (
                """(HLT_Ele35_WPTight_Gsf == true)
                && (first_iso_electron_pt > 38)
                && (second_iso_electron_pt > 15)"""
        )
        events = events.Define(
            "trig_emu_MET",
            (
                """subtrigger_MET
                && ((first_iso_electron_pt > 25) && (first_iso_muon_pt > 15))
                || ((first_iso_electron_pt > 15) && (first_iso_muon_pt > 25))"""
        )
        events = events.Define(
            "trig_mumu_MET",
            """(subtrigger_MET && (first_iso_muon_pt > 25) && (second_iso_muon_pt > 15))"""
        )
        events = events.Define(
            "trig_ee_MET", 
                """subtrigger_MET
                && (first_iso_electron_pt > 25)
                && (second_iso_electron_pt > 15)"""
        )
        if is_mc:
            events = events.Define(
                "trig_emu", "(trig_emu_MuonEG1 || trig_emu_MuonEG2 || trig_emu_SingleMuon || trig_emu_SingleElectron)"
            )
            events = events.Define(
                "trig_mumu", "(trig_mumu_DoubleMuon || trig_mumu_SingleMuon)"
            )
            events = events.Define(
                "trig_ee", "(trig_ee_DoubleEG || trig_ee_SingleElectron"
            )
            # TODO:
            # Insert dilepton events triggered on MET, NO overlap removal in data - use these for trigger scale factors
            # So need events passing dilepton cuts + MET trigger/MET cut as denominator of SFs,
            # and events passing dilepton cuts + MET trigger/MET cut + dilepton triggers for numerator -> Get efficiency
        else:
            if data_stream == "MuonEG":
                # Select emu if they pass any triggers, we'll pick up the rest from the exclusive events in SingleMuon and SingleElectron
                events = events.Define(
                    "trig_emu", "trig_emu_MuonEG1 || trig_emu_MuonEG2"
                )# , "trig_emu_SingleMuon", "trig_emu_SingleElectron"), #This doesn't veto properly
                events = events.Define("trig_mumu", "false")
                events = events.Define("trig_ee", "false")
            elif data_stream == "DoubleMuon":
                events = events.Define("trig_emu", "false")
                events = events.Define(
                    "trig_mumu", "trig_mumu_DoubleMuon"
                )# , "trig_mumu_SingleMuon"),
                events = events.Define("trig_ee", "false")
            elif data_stream == "SingleMuon":
                events = events.Define(
                    "trig_emu", "(trig_emu_MuonEG1 == false && trig_emu_MuonEG2 == false && trig_emu_SingleMuon == true)"
                )
                events = events.Define(
                    "trig_mumu", "(trig_mumu_DoubleMuon == false && trig_mumu_SingleMuon == true)"
                )
                events = events.Define("trig_ee", "false")
            elif data_stream == "DoubleEG":
                events = events.Define("trig_emu", "false")
                events = events.Define("trig_mumu", "false")
                events = events.Define(
                    "trig_ee", "trig_ee_DoubleEG"
                )# , "trig_ee_SingleElectron"),
            elif data_stream == "SingleElectron":
                events = events.Define(
                    "trig_emu", 
                    """trig_emu_MuonEG1 == false && trig_emu_MuonEG2 == false && 
                    trig_emu_SingleMuon == false && trig_emu_SingleElectron == true)"""
                )
                events = events.Define("trig_mumu", "false")
                events = events.Define(
                    "trig_ee", "(trig_ee_DoubleEG == false && trig_ee_SingleElectron == true)"
                )
            elif data_stream == "MET":
                raise NotImplementedError
            else:
                raise ValueError(
                    f"data_stream must be specified for is_mc=False (for event overlap removal), input: {data_stream}"
                )
            # Do e.g. trig_masks.require(trig_emu_MuonEG1=False, trig_emu_MuonEG2=False, trig_emu_SingleMuon=True)
            # for the SingleMuon datastream/backup trigger, so no double counting events in the MuonEG datastream+trigger
    elif era == "2016":
        raise NotImplementedError
    return trig_masks
