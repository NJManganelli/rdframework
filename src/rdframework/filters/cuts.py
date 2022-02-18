from __future__ import annotations

from typing import Any


def PV_MET_filter(
    events: Any,
    era: str,
    is_mc: bool,
    min_NDoF: int = 4,
    max_abs_z: float = 24.0,
    max_rho: float = 2.0,
    is_ultra_legacy: bool = True,
    is_fastsim_MC: bool = False,
    include_HF: bool = False,
    return_applied_flags: bool = False,
) -> Any:
    if is_fastsim_MC is True:
        raise NotImplementedError(
            "Appropriate settings for FastSim MonteCarlo is not yet included in this function"
        )
    flags = []
    if is_ultra_legacy is True:
        # Use the full Flag_* name just to match better with the twiki, strip off the prefix later, not too expensive
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#UL_data as of Nov 15, 2021
        if era in ["2016", "2017", "2018"]:
            flags.append("Flag_goodVertices")
            flags.append("Flag_globalSuperTightHalo2016Filter")
            flags.append("Flag_HBHENoiseFilter")
            flags.append("Flag_HBHENoiseIsoFilter")
            flags.append("Flag_EcalDeadCellTriggerPrimitiveFilter")
            flags.append("Flag_BadPFMuonFilter")
            flags.append("Flag_BadPFMuonDzFilter")
            flags.append("Flag_eeBadScFilter")
        if era in ["2017", "2018"]:
            flags.append("Flag_ecalBadCalibFilter")
            if include_HF is True:
                flags.append("Flag_hfNoisyHitsFilter")
        # Not recommended:
        # "Flag_BadChargedCandidateFilter"
    elif is_ultra_legacy is False:
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2018_data_EOY as of Nov 15, 2021
        if era in ["2016", "2017", "2018"]:
            flags.append("Flag_goodVertices")
            flags.append("Flag_globalSuperTightHalo2016Filter")
            flags.append("Flag_HBHENoiseFilter")
            flags.append("Flag_HBHENoiseIsoFilter")
            flags.append("Flag_EcalDeadCellTriggerPrimitiveFilter")
            flags.append("Flag_BadPFMuonFilter")
            if is_mc is False:
                flags.append("Flag_eeBadScFilter")
        if era in ["2017", "2018"]:
            flags.append("Flag_ecalBadCalibFilterV2")
    else:
        raise NotImplementedError(
            "legacy should be True or False, True if it's an (Ultra-)Legacy sample."
        )

    PV_mask = f"return (PV_ndof >= {min_NDoF}) && (abs(PV_z) < {max_abs_z}) && (sqrt(PV_x * PV_x + PV_y * PV_y) < {max_rho});"
    PV_nicename = "PV Filters: " + PV_mask.replace("&&", "and").replace("PV_", "")

    flag_mask = "return "
    flag_nicename = "MET Filters: "
    for i, raw_mask in enumerate(flags):
        if 0 < i < len(flags):
            flag_mask += " && "
        flag_mask += raw_mask
        flag_nicename += raw_mask[5:]
    flag_mask += ";"
    filtered = events.Filter(PV_mask, PV_nicename).Filter(flag_mask, flag_nicename)
    if return_applied_flags:
        return filtered, flags
    else:
        return filtered
