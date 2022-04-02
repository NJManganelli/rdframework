from __future__ import annotations

from pathlib import Path
from typing import Any

import ROOT

if not hasattr(ROOT, "METXYCorr_Met_MetPhi"):
    src = Path(__file__).parent / "met.cpp"
    ROOT.gROOT.ProcessLine(f".L {src}")


def MET_xy_corrector(
    input_df: Any,
    era: str,
    is_mc: bool,
    MET_fields: list[str] | None = None,
    MET_phi_fields: list[str] | None = None,
    npv_fields: list[str] | None = None,
    run_period: str | None = None,
    is_ultra_legacy: bool = True,
    pre_post_VFP: str | None = None,
    is_fastsim_MC: bool = False,
    is_PUPPI: bool = False,
    # run_branch = "run",
) -> Any:

    # Set defaults
    input_MET_fields = MET_fields if MET_fields is not None else ["MET", "pt"]
    input_MET_phi_fields = (
        MET_phi_fields if MET_phi_fields is not None else ["MET", "phi"]
    )
    input_npv_fields = npv_fields if npv_fields is not None else ["PV", "npvs"]

    rdf = input_df

    uncormet = "_".join(input_MET_fields)
    cormet = "MET_rdf_xycorr_pt"
    uncormet_phi = "_".join(input_MET_phi_fields)
    cormet_phi = "MET_rdf_xycorr_phi"
    runnb = "run"
    isMC = str(is_mc).lower()
    npv = "_".join(input_npv_fields)
    isUL = str(is_ultra_legacy).lower()
    ispuppi = str(is_PUPPI).lower()
    year = str(era)
    if is_ultra_legacy and era == "2016":
        if pre_post_VFP == "preVFP":
            year += "APV"
        elif pre_post_VFP == "postVFP":
            year += "nonAPV"
        else:
            raise ValueError(
                f"Invalid choice of pre_post_VFP ({pre_post_VFP}) for year ({era})."
            )

    if is_mc:
        call = f'return METXYCorr_Met_MetPhi({uncormet}, {uncormet_phi}, {runnb}, "{year}", {isMC}, {npv}+1, {isUL}, {ispuppi});'
    else:
        call = f'return METXYCorr_Met_MetPhi({uncormet}, {uncormet_phi}, {runnb}, "{year}", {isMC}, {npv}, {isUL}, {ispuppi});'
    # METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, TString year, bool isMC, int npv, bool isUL =false,bool ispuppi=false)

    met_doublet = "MET_xycorr_doublet"
    rdf = rdf.Define(met_doublet, call)
    rdf = rdf.Define(cormet, f"{met_doublet}.first")
    rdf = rdf.Define(cormet_phi, f"{met_doublet}.second")
    return rdf
