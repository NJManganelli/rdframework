from __future__ import annotations

from typing import Any


def vidUnpackedWP(electrons: Any) -> Any:
    """Return a dictionary of the cuts in the electron cutBasedID,
    e.g. results["GsfEleEInverseMinusPInverseCut"] will be 0 (fail), 1, 2, 3, or 4 (tight)"""
    results = dict()
    for name, shift in zip(
        [
            "MinPtCut",
            "GsfEleSCEtaMultiRangeCut",
            "GsfEleDEtaInSeedCut",
            "GsfEleDPhiInCut",
            "GsfEleFull5x5SigmaIEtaIEtaCut",
            "GsfEleHadronicOverEMEnergyScaledCut",
            "GsfEleEInverseMinusPInverseCut",
            "GsfEleRelPFIsoScaledCut",
            "GsfEleConversionVetoCut",
            "GsfEleMissingHitsCut",
        ],
        range(0, 28, 3),
    ):
        results[name] = (electrons.vidNestedWPBitmap >> shift) & 0b111
    return results


def vidUnpackedWPSelection(electrons: Any, level: int) -> Any:
    """Return a dictionary of boolean masks for the electron cutBasedID,
    e.g. results["GsfEleEInverseMinusPInverseCut"] will be True if the result value is >= level"""
    results = dict()
    for name, cut_level in vidUnpackedWP(electrons).items():
        results[name] = cut_level >= level
    return results


def select_electrons_cutBased(
    electrons: Any,
    collection_name: str,
    columns: list[str] | None,
    electron_id: int | str,
    min_pt: float,
    max_eta: float = 2.5,
    max_ip3d_barrel: float = 0.05,
    max_ip3d_endcap: float = 0.10,
    max_dz_barrel: float = 0.1,
    max_dz_endcap: float = 0.2,
    return_mask: bool = True,
    invert_cuts: list[str] | None = None,
) -> Any:
    if isinstance(electron_id, str):
        if electron_id.lower() == "fail":
            e_id = 0
        elif electron_id.lower() == "veto":
            e_id = 1
        elif electron_id.lower() == "loose":
            e_id = 2
        elif electron_id.lower() == "medium":
            e_id = 3
        elif electron_id.lower() == "tight":
            e_id = 4
        else:
            raise ValueError(f"{electron_id} is not a supported cutBased electron ID")
    else:
        e_id = electron_id
    # We apply the EGamma recommendations
    mask = (
        (np.abs(electrons.eta) < 1.4442)
        & (np.abs(electrons.ip3d) < max_ip3d_barrel)
        & (np.abs(electrons.dz) < max_dz_barrel)
    ) | (
        (np.abs(electrons.eta) > 1.5660)
        & (np.abs(electrons.eta) <= max_eta)
        & (np.abs(electrons.ip3d) < max_ip3d_endcap)
        & (np.abs(electrons.dz) < max_dz_endcap)
    )
    mask = mask & (electrons.pt >= min_pt)
    if isinstance(invert_cuts, list) and len(invert_cuts) > 0:
        for name, pass_fail in vidUnpackedWPSelection(electrons, level=e_id).items():
            if name not in invert_cuts:
                # Need the cut to have passed at the minimum level
                mask = mask & pass_fail
            else:
                # Want the cut to have failed, so ~pass_fail will be True
                mask = mask & (~pass_fail)
    else:
        mask = mask & (electrons.cutBased >= e_id)

    if return_mask == True:
        return electrons[mask], mask
    else:
        return electrons[mask]


def select_muons_cutBased(
    muons: Any,
    collection_name: str,
    columns: list[str] | None,
    muon_id: str,
    muon_iso: int | str,
    min_pt: float,
    max_eta: float = 2.4,
    max_ip3d: float = 0.10,
    max_dz: float = 0.2,
    return_mask: bool = True,
    invert_iso: bool = False,
) -> Any:
    if isinstance(muon_iso, str):
        # 1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight
        if muon_iso.lower() in ["fail", "pfisoveryveryloose"]:
            m_iso = 0
        elif muon_iso.lower() in ["pfisoveryloose", "veryloose"]:
            m_iso = 1
        elif muon_iso.lower() in ["pfisoloose", "loose"]:
            m_iso = 2
        elif muon_iso.lower() in ["pfisomedium", "medium"]:
            m_iso = 3
        elif muon_iso.lower() in ["pfisotight", "tight"]:
            m_iso = 4
        elif muon_iso.lower() in ["pfisoverytight", "verytight"]:
            m_iso = 5
        elif muon_iso.lower() in ["pfisoveryverytight", "veryverytight"]:
            m_iso = 6
        else:
            raise ValueError(f"{muon_iso} is not a supported cutBased muon ID")
    else:
        m_iso = muon_iso

    # We apply the cuts
    mask = (
        (muons.pt >= min_pt)
        & (np.abs(muons.eta) <= max_eta)
        & (np.abs(muons.ip3d) < max_ip3d)
        & (np.abs(muons.dz) < max_dz)
    )
    if invert_iso == True:
        mask = mask & (muons.pfIsoId < m_iso)
    else:
        mask = mask & (muons.pfIsoId >= m_iso)

    # apply the cutbased ID, which is stored in separate boolean branches/fields
    if muon_id.lower() == "fail":
        pass  # accept all the muons...
    elif muon_id.lower() == "veto":
        raise ValueError("Muons do not have a veto cutBased working point")
    elif muon_id.lower() == "loose":
        mask = mask & muons.looseId
    elif muon_id.lower() == "medium":
        mask = mask & muons.mediumId
    elif muon_id.lower() == "tight":
        mask = mask & muons.tightId

    if return_mask == True:
        return muons[mask], mask
    else:
        return muons[mask]
