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
    if not noniso_muons.endswith("_"):
        noniso_muons += "_"
    if not noniso_electrons.endswith("_"):
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
    iso_muons: Any,
    iso_electrons: Any,
    era: str,
    is_mc: bool,
    run_period: str | None = None,
    data_stream: str | None = None,
    noniso_muons: Any | None = None,
    noniso_electrons: Any | None = None,
) -> Any:
    # This needs to be rethought. Maybe go down to the TriggerObject level and rely upon that instead!
    # Leading lepton pt's for trigger matching, broadcast-able with trigger arrays (1D in events axis)
    first_iso_muon_pt = ak.flatten(
        select_with_default(iso_muons, "pt", elements=0, default=0.0, axis=-1)
    )
    second_iso_muon_pt = ak.flatten(
        select_with_default(iso_muons, "pt", elements=1, default=0.0, axis=-1)
    )
    first_iso_electron_pt = ak.flatten(
        select_with_default(iso_electrons, "pt", elements=0, default=0.0, axis=-1)
    )
    second_iso_electron_pt = ak.flatten(
        select_with_default(iso_electrons, "pt", elements=1, default=0.0, axis=-1)
    )

    # a_first_iso_muon_pt = ak.flatten(ak.fill_none(ak.pad_none(iso_muons[:, 0:1].pt, target=1, axis=-1), value=0.0, axis=-1))
    # a_second_iso_muon_pt = ak.flatten(ak.fill_none(ak.pad_none(iso_muons[:, 1:2].pt, target=1, axis=-1), value=0.0, axis=-1))
    # a_first_iso_electron_pt = ak.flatten(ak.fill_none(ak.pad_none(iso_electrons[:, 0:1].pt, target=1, axis=-1), value=0.0, axis=-1))
    # a_second_iso_electron_pt = ak.flatten(ak.fill_none(ak.pad_none(iso_electrons[:, 1:2].pt, target=1, axis=-1), value=0.0, axis=-1))

    # assert ak.all(first_iso_muon_pt == a_first_iso_muon_pt)
    # assert ak.all(second_iso_muon_pt == a_second_iso_muon_pt)
    # assert ak.all(first_iso_electron_pt == a_first_iso_electron_pt)
    # assert ak.all(second_iso_electron_pt == a_second_iso_electron_pt)

    trig_masks = PackedSelection(dtype="uint32")
    # Store masks as trig_channel_datastream<number>, e.g. trig_emu_MuonEG1 and ...MuonEG2 are the two electron-muon
    # dilepton trigger paths, selecting events from the MuonEG datastream if is_mc=False and being used as a veto from
    # other datastreams... in 2018, the DoubleElectron datastream doesn't exist anymore, so... care must be taken
    # minimally dilepton events... does not include the single-lepton events, with or without a noniso 2nd lepton
    # WAITWAITWAIT... does it still make sense to make these 'dilepton' triggers? hard to say with ML events... "\
    # "yeah, probably, then make another one for single lepton or ML events, and require appropriate lepton pts")
    # Is there a good reason why everyone doesn't use the real trigger-level objects? That would be cleaner...
    subtrigger_MET = (
        events.HLT.PFMETTypeOne200_HBHE_BeamHaloCleaned
        | events.HLT.PFMET200_HBHECleaned
        | events.HLT.PFMET200_NotCleaned
    ) & (events.MET.pt > 210)
    false_val = (
        events.MET.pt < 0
    )  # always false, for filling non-permissible data_stream + trigger combos
    if era == "2018":
        trig_masks.add(
            "trig_emu_MuonEG1",
            (
                (events.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == True)
                & (first_iso_electron_pt > 25)
                & (first_iso_muon_pt > 15)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_emu_MuonEG2",
            (
                (events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == True)
                & (first_iso_electron_pt > 15)
                & (first_iso_muon_pt > 25)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_mumu_DoubleMuon",
            (
                (events.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 == True)
                & (first_iso_muon_pt > 25)
                & (second_iso_muon_pt > 15)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_ee_EGamma1",
            (
                (events.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL == True)
                & (first_iso_electron_pt > 25)
                & (second_iso_electron_pt > 15)
            ),
            fill_value=False,
        )
        # backup channels not used in the current analysis, but gains back some events...
        trig_masks.add(
            "trig_emu_SingleMuon",
            (
                (events.HLT.IsoMu24 == True)
                & (first_iso_electron_pt > 15)
                & (first_iso_muon_pt > 27)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_emu_EGamma",
            (
                (events.HLT.Ele32_WPTight_Gsf == True)
                & (first_iso_electron_pt > 35)
                & (first_iso_muon_pt > 15)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_mumu_SingleMuon",
            (
                (events.HLT.IsoMu24 == True)
                & (first_iso_muon_pt > 27)
                & (second_iso_muon_pt > 15)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_ee_EGamma2",
            (
                (events.HLT.Ele32_WPTight_Gsf == True)
                & (first_iso_electron_pt > 35)
                & (second_iso_electron_pt > 15)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_emu_MET",
            (
                subtrigger_MET
                & ((first_iso_electron_pt > 25) & (first_iso_muon_pt > 15))
                | ((first_iso_electron_pt > 15) & (first_iso_muon_pt > 25))
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_mumu_MET",
            (subtrigger_MET & (first_iso_muon_pt > 25) & (second_iso_muon_pt > 15)),
            fill_value=False,
        )
        trig_masks.add(
            "trig_ee_MET",
            (
                subtrigger_MET
                & (first_iso_electron_pt > 25)
                & (second_iso_electron_pt > 15)
            ),
            fill_value=False,
        )
        if is_mc:
            trig_masks.add(
                "trig_emu",
                trig_masks.any(
                    "trig_emu_MuonEG1",
                    "trig_emu_MuonEG2",
                    "trig_emu_SingleMuon",
                    "trig_emu_EGamma",
                ),
                fill_value=False,
            )
            trig_masks.add(
                "trig_mumu",
                trig_masks.any("trig_mumu_DoubleMuon", "trig_mumu_SingleMuon"),
                fill_value=False,
            )
            trig_masks.add(
                "trig_ee",
                trig_masks.any("trig_ee_EGamma1", "trig_ee_EGamma2"),
                fill_value=False,
            )
        else:
            if data_stream == "MuonEG":
                # Select emu if they pass any triggers, we'll pick up the rest from the exclusive events in SingleMuon and EGamma
                trig_masks.add(
                    "trig_emu",
                    trig_masks.any(
                        "trig_emu_MuonEG1", "trig_emu_MuonEG2"
                    ),  # , "trig_emu_SingleMuon", "trig_emu_EGamma"),
                    fill_value=False,
                )
                trig_masks.add("trig_mumu", false_val, fill_value=False)
                trig_masks.add("trig_ee", false_val, fill_value=False)
            elif data_stream == "DoubleMuon":
                trig_masks.add("trig_emu", false_val, fill_value=False)
                trig_masks.add(
                    "trig_mumu",
                    trig_masks.any(
                        "trig_mumu_DoubleMuon"
                    ),  # , "trig_mumu_SingleMuon"),
                    fill_value=False,
                )
                trig_masks.add("trig_ee", false_val, fill_value=False)
            elif data_stream == "SingleMuon":
                trig_masks.add(
                    "trig_emu",
                    trig_masks.require(
                        trig_emu_MuonEG1=False,
                        trig_emu_MuonEG2=False,
                        trig_emu_SingleMuon=True,
                    ),
                    fill_value=False,
                )
                trig_masks.add(
                    "trig_mumu",
                    trig_masks.require(
                        trig_mumu_DoubleMuon=False, trig_mumu_SingleMuon=True
                    ),
                    fill_value=False,
                )
                trig_masks.add("trig_ee", false_val, fill_value=False)
            elif data_stream == "EGamma":
                trig_masks.add(
                    "trig_emu",
                    trig_masks.require(
                        trig_emu_MuonEG1=False,
                        trig_emu_MuonEG2=False,
                        trig_emu_SingleMuon=False,
                        trig_emu_EGamma=True,
                    ),
                    fill_value=False,
                )
                trig_masks.add("trig_mumu", false_val, fill_value=False)
                trig_masks.add(
                    "trig_ee",
                    trig_masks.any("trig_ee_EGamma1", "trig_ee_EGamma2"),
                    fill_value=False,
                )
            elif data_stream == "MET":
                raise NotImplementedError
            else:
                raise ValueError(
                    f"data_stream must be specified for is_mc=False (for event overlap removal), input: {data_stream}"
                )
            # Do e.g. trig_masks.require(trig_emu_MuonEG1=False, trig_emu_MuonEG2=False, trig_emu_SingleMuon=True)
            # for the SingleMuon datastream/backup trigger, so no double counting events in the MuonEG datastream+trigger
    elif era == "2017":
        trig_masks.add(
            "trig_emu_MuonEG1",
            (
                (events.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == True)
                & (first_iso_electron_pt > 25)
                & (first_iso_muon_pt > 15)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_emu_MuonEG2",
            (
                (events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == True)
                & (first_iso_electron_pt > 15)
                & (first_iso_muon_pt > 25)
            ),
            fill_value=False,
        )
        if not is_mc and run_period == "B":
            trig_masks.add(
                "trig_mumu_DoubleMuon",
                (
                    (events.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ == True)
                    & (first_iso_muon_pt > 25)
                    & (second_iso_muon_pt > 15)
                ),
                fill_value=False,
            )
        else:
            trig_masks.add(
                "trig_mumu_DoubleMuon",
                (
                    (events.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 == True)
                    & (first_iso_muon_pt > 25)
                    & (second_iso_muon_pt > 15)
                ),
                fill_value=False,
            )
        trig_masks.add(
            "trig_ee_DoubleEG",
            (
                (events.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL == True)
                & (first_iso_electron_pt > 25)
                & (second_iso_electron_pt > 15)
            ),
            fill_value=False,
        )

        trig_masks.add(
            "trig_emu_SingleMuon",
            (
                (events.HLT.IsoMu27 == True)
                & (first_iso_electron_pt > 15)
                & (first_iso_muon_pt > 30)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_emu_SingleElectron",
            (
                (events.HLT.Ele35_WPTight_Gsf == True)
                & (first_iso_electron_pt > 38)
                & (first_iso_muon_pt > 15)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_mumu_SingleMuon",
            (
                (events.HLT.IsoMu27 == True)
                & (first_iso_muon_pt > 30)
                & (second_iso_muon_pt > 15)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_ee_SingleElectron",
            (
                (events.HLT.Ele35_WPTight_Gsf == True)
                & (first_iso_electron_pt > 38)
                & (second_iso_electron_pt > 15)
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_emu_MET",
            (
                subtrigger_MET
                & ((first_iso_electron_pt > 25) & (first_iso_muon_pt > 15))
                | ((first_iso_electron_pt > 15) & (first_iso_muon_pt > 25))
            ),
            fill_value=False,
        )
        trig_masks.add(
            "trig_mumu_MET",
            (subtrigger_MET & (first_iso_muon_pt > 25) & (second_iso_muon_pt > 15)),
            fill_value=False,
        )
        trig_masks.add(
            "trig_ee_MET",
            (
                subtrigger_MET
                & (first_iso_electron_pt > 25)
                & (second_iso_electron_pt > 15)
            ),
            fill_value=False,
        )
        if is_mc:
            trig_masks.add(
                "trig_emu",
                trig_masks.any(
                    "trig_emu_MuonEG1",
                    "trig_emu_MuonEG2",
                    "trig_emu_SingleMuon",
                    "trig_emu_SingleElectron",
                ),
                fill_value=False,
            )
            trig_masks.add(
                "trig_mumu",
                trig_masks.any("trig_mumu_DoubleMuon", "trig_mumu_SingleMuon"),
                fill_value=False,
            )
            trig_masks.add(
                "trig_ee",
                trig_masks.any("trig_ee_DoubleEG", "trig_ee_SingleElectron"),
                fill_value=False,
            )
            # TODO:
            # Insert dilepton events triggered on MET, NO overlap removal in data - use these for trigger scale factors
            # So need events passing dilepton cuts + MET trigger/MET cut as denominator of SFs,
            # and events passing dilepton cuts + MET trigger/MET cut + dilepton triggers for numerator -> Get efficiency
        else:
            if data_stream == "MuonEG":
                # Select emu if they pass any triggers, we'll pick up the rest from the exclusive events in SingleMuon and SingleElectron
                trig_masks.add(
                    "trig_emu",
                    trig_masks.any(
                        "trig_emu_MuonEG1", "trig_emu_MuonEG2"
                    ),  # , "trig_emu_SingleMuon", "trig_emu_SingleElectron"), #This doesn't veto properly
                    fill_value=False,
                )
                trig_masks.add("trig_mumu", false_val, fill_value=False)
                trig_masks.add("trig_ee", false_val, fill_value=False)
            elif data_stream == "DoubleMuon":
                trig_masks.add("trig_emu", false_val, fill_value=False)
                trig_masks.add(
                    "trig_mumu",
                    trig_masks.any(
                        "trig_mumu_DoubleMuon"
                    ),  # , "trig_mumu_SingleMuon"),
                    fill_value=False,
                )
                trig_masks.add("trig_ee", false_val, fill_value=False)
            elif data_stream == "SingleMuon":
                trig_masks.add(
                    "trig_emu",
                    trig_masks.require(
                        trig_emu_MuonEG1=False,
                        trig_emu_MuonEG2=False,
                        trig_emu_SingleMuon=True,
                    ),
                    fill_value=False,
                )
                trig_masks.add(
                    "trig_mumu",
                    trig_masks.require(
                        trig_mumu_DoubleMuon=False, trig_mumu_SingleMuon=True
                    ),
                    fill_value=False,
                )
                trig_masks.add("trig_ee", false_val, fill_value=False)
            elif data_stream == "DoubleEG":
                trig_masks.add("trig_emu", false_val, fill_value=False)
                trig_masks.add("trig_mumu", false_val, fill_value=False)
                trig_masks.add(
                    "trig_ee",
                    trig_masks.any("trig_ee_DoubleEG"),  # , "trig_ee_SingleElectron"),
                    fill_value=False,
                )
            elif data_stream == "SingleElectron":
                trig_masks.add(
                    "trig_emu",
                    trig_masks.require(
                        trig_emu_MuonEG1=False,
                        trig_emu_MuonEG2=False,
                        trig_emu_SingleMuon=False,
                        trig_emu_SingleElectron=True,
                    ),
                    fill_value=False,
                )
                trig_masks.add("trig_mumu", false_val, fill_value=False)
                trig_masks.add(
                    "trig_ee",
                    trig_masks.require(
                        trig_ee_DoubleEG=False, trig_ee_SingleElectron=True
                    ),
                    fill_value=False,
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
