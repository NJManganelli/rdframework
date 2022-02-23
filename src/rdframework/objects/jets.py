from __future__ import annotations

from typing import Any


def btag_jets(
    events: Any,
    btagger: str,
    WP: str,
    era: str,
    is_ultra_legacy: bool = True,
    pre_post_VFP: str | None = None,
    return_btag_dict: bool = False,
) -> Any:
    # Nested btag WP dictionary of the form dict[is_ultra_legacy][era][btagger]
    # From https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation#Recommendation_for_13_TeV_Data
    bTagWorkingPointDict = dict()
    bTagWorkingPointDict[False] = {
        "2016": {
            "DeepCSV": {"L": 0.2217, "M": 0.6321, "T": 0.8953, "Var": "btagDeepB"},
            "DeepJet": {"L": 0.0614, "M": 0.3093, "T": 0.7221, "Var": "btagDeepFlavB"},
        },
        "2017": {
            "CSVv2": {"L": 0.5803, "M": 0.8838, "T": 0.9693, "Var": "btagCSVV2"},
            "DeepCSV": {"L": 0.1522, "M": 0.4941, "T": 0.8001, "Var": "btagDeepB"},
            "DeepJet": {"L": 0.0521, "M": 0.3033, "T": 0.7489, "Var": "btagDeepFlavB"},
        },
        "2018": {
            "DeepCSV": {"L": 0.1241, "M": 0.4184, "T": 0.7527, "Var": "btagDeepB"},
            "DeepJet": {"L": 0.0494, "M": 0.2770, "T": 0.7264, "Var": "btagDeepFlavB"},
        },  # Non-UL WPs
    }
    bTagWorkingPointDict[True] = {
        "2016preVFP": {
            "DeepCSV": {"L": 0.2027, "M": 0.6001, "T": 0.8819, "Var": "btagDeepB"},
            "DeepJet": {"L": 0.0508, "M": 0.2598, "T": 0.6502, "Var": "btagDeepFlavB"},
        },
        "2016postVFP": {
            "DeepCSV": {"L": 0.1918, "M": 0.5847, "T": 0.8767, "Var": "btagDeepB"},
            "DeepJet": {"L": 0.0480, "M": 0.2489, "T": 0.6377, "Var": "btagDeepFlavB"},
        },
        "2017": {
            "DeepCSV": {"L": 0.1355, "M": 0.4506, "T": 0.7738, "Var": "btagDeepB"},
            "DeepJet": {"L": 0.0532, "M": 0.3040, "T": 0.7476, "Var": "btagDeepFlavB"},
        },
        "2018": {
            "DeepCSV": {"L": 0.1208, "M": 0.4168, "T": 0.7665, "Var": "btagDeepB"},
            "DeepJet": {"L": 0.0490, "M": 0.2783, "T": 0.7100, "Var": "btagDeepFlavB"},
        },  # UL WPs
    }

    # Better to do 2016 -> 2016preVFP 2016postVFP?
    if is_ultra_legacy and era == "2016":
        if pre_post_VFP in ["preVFP", "postVFP"]:
            key_era = era + pre_post_VFP
        else:
            raise ValueError("UltraLegacy 2016 requires a 'preVFP' or 'postVFP' tag")
    else:
        key_era = era
    subset = bTagWorkingPointDict[is_ultra_legacy][key_era][btagger]
    btagVar, btagWP = str(subset["Var"]), subset[WP]
    results = [events.Define("btagmask", f"return Jet_{btagVar} >= {btagWP};")]
    if return_btag_dict:
        results.append(subset)
    if len(results) > 1:
        return results
    else:
        return results[0]


def select_jets(
    events: Any,
    input_collection: str,
    output_collection: str,
    columns: list[str] | str | None,
    isolated_leptons: list[str] | None,
    clean_algo_or_dR: float | str,
    jet_min_pt: float,
    jet_max_eta: float,
    jet_id: str,
    jet_pu_id: str | None = None,
    btagging_configuration: dict[str, Any] | None = None,
    era: str | None = None,
    is_ultra_legacy: bool = True,
    pre_post_VFP: str | None = None,
    fix_inverted_pu_id_bits: bool = False,
    sort_column: str | None = "pt",
    sort_reverse: bool = True,
) -> Any:
    """pass in RDataFrame and list of names of isolated leptons to clean against. can use PFMatching or deltaR

    General notes: for tracker acceptance, max_eta of 2.5 for 2017-2018, 2.4 for 2016
    Jet PU ID is good to apply for jets with 20 GeV < pt < 50 GeV, especially outside the tracker acceptance
    Due to a configuration mistake, one or both of the UltraLegacy 2016 samples have the Jet PU ID bits reversed!
    Also recommended to apply PU ID Loose if using the DeepJet tagger, as it was trained with that WP in place
    In 2016, 'loose" jet ID is available, but 2017 and 2018 only have 'tight' and 'tightlepveto'

    """
    # This function either needs to be aware of jet JES/JER variations in pt, or will be called multiple times!
    # The latter might be extremely costly and inefficient if you have to use e.g. all ~20 JEC variations!
    # Once the Vary function is around, that should be called on the pt, then the variations should work downstream for
    # anything using jet pt
    if not input_collection.endswith("_"):
        input_collection += "_"
    if not output_collection.endswith("_"):
        output_collection += "_"

    if fix_inverted_pu_id_bits:
        raise NotImplementedError(
            "Patching of the inverted Jet PU ID bits in 2016 UL not implemented"
        )

    if jet_id.lower() in ["loose", "l"]:
        jet_min_id = 1  # don't use in 2017/2018!
    elif jet_id.lower() in ["tight", "t"]:
        jet_min_id = 2
    elif jet_id.lower() in ["tightlepveto", "tlv"]:
        jet_min_id = 6
    else:
        raise ValueError(f"Unsupported jet_id value{jet_id}")

    # This presumes _pt is a Vary'd column, all systematics accounted for!
    mask = (
        f"auto jmask = ({input_collection}pt > {jet_min_pt}) && abs({input_collection}eta) <= {jet_max_eta}"
        f"&& ({input_collection}jetId >= {jet_min_id})"
    )

    if jet_pu_id:
        if jet_pu_id.lower() in ["loose", "l"]:
            jet_min_pu_id = 4
        elif jet_pu_id.lower() in ["medium", "m"]:
            jet_min_pu_id = 6
        elif jet_pu_id.lower() in ["tight", "t"]:
            jet_min_pu_id = 7
            # jet_mask[syst_name] = jet_mask[syst_name] & ( (getattr(jets, "pt_" + syst_variation) > 50.0) | jets.puId == 7)
        else:
            raise ValueError("Invalid Jet PU Id selected")
        mask = (
            mask
            + f" && (({input_collection}pt > 50.0) || ({input_collection}puId >= {jet_min_pu_id}))"
        )
    mask += ";\n"

    # allow for 0 to many isolated lepton collections to be passed in...
    if isolated_leptons:
        for lep_collection in isolated_leptons:
            if isinstance(clean_algo_or_dR, str):  # PFMatching
                # FIXME: add this to logging with a check against DefinedColumnNames
                # WARNING: PFMatching against a reduced collection will produce incorrect cross-cleaning
                mask = mask + (
                    f"for(int i=0; i < {lep_collection}_jetIdx.size(); ++i){{\n"
                    f"  jmask = jmask && ({input_collection}idx != {lep_collection}_jetIdx.at(i));\n"
                    "}\n"
                )
            elif isinstance(clean_algo_or_dR, float):  # DeltaR
                mask = mask + (
                    f"for(int i=0; i < {lep_collection}_jetIdx.size(); ++i){{\n"
                    "  ROOT::VecOps::RVec<double> dr;\n"
                    "  for(int j=0; j < jmask.size(); ++j){\n"
                    f"    dr.push_back(ROOT::VecOps::DeltaR({input_collection}eta.at(j),\n"
                    f"                                      {lep_collection}_eta.at(i),\n"
                    f"                                      {input_collection}phi.at(j),\n"
                    f"                                      {lep_collection}_phi.at(i)\n"
                    "                       )\n"
                    "                );\n"
                    "  } // end loop over jmask\n"
                    f"  jmask = jmask && dr >= {clean_algo_or_dR};\n"
                    "  dr.clear();\n"
                    "} //end loop over lep_collection\n"
                )
    mask = mask + "return jmask;"

    avail_columns = [
        str(col)
        for col in events.GetColumnNames()
        if str(col).startswith(input_collection)
    ]
    if f"{input_collection}idx" not in avail_columns:
        events = events.Define(
            f"{input_collection}idx",
            f"ROOT::VecOps::RVec<UInt_t> empty = {{}}; return {avail_columns[0]}.size() > 0 ? ROOT::VecOps::RVec<UInt_t>(Combinations({avail_columns[0]}, 1).at(0)) : empty;",
        )
    events = events.Define(f"{output_collection}jetmask", mask)
    if btagging_configuration is not None:
        btagger = str(btagging_configuration.get("btagger"))
        WP = str(btagging_configuration.get("WP"))
        era = str(btagging_configuration.get("era", era))
        is_ultra_legacy = btagging_configuration.get("is_ultra_legacy", is_ultra_legacy)
        pre_post_VFP = btagging_configuration.get("is_ultra_legacy", pre_post_VFP)
        events = btag_jets(
            events,
            btagger=btagger,
            WP=WP,
            era=era,
            is_ultra_legacy=is_ultra_legacy,
            pre_post_VFP=pre_post_VFP,
            return_btag_dict=False,
        )
    if isinstance(columns, list):
        sel_columns = [
            col[len(input_collection) :]
            for col in columns
            if col.startswith(f"{input_collection}")
        ]
    elif isinstance(columns, str):
        raise NotImplementedError("regexp not currently supported")
    else:
        sel_columns = [col[len(input_collection) :] for col in avail_columns] + [
            f"{input_collection}idx"
        ]

    if sort_column:
        if sort_reverse:
            events = events.Define(
                f"{output_collection}jettake",
                f"return Reverse(Argsort({input_collection}{sort_column}[{output_collection}jetmask]));",
            )
        else:
            events = events.Define(
                "jettake",
                f"return Argsort({input_collection}{sort_column}[{output_collection}jetmask]);",
            )
    else:
        events = events.Define(
            f"{output_collection}jettake",
            f"return {input_collection}idx[{output_collection}jetmask];",
        )
    events = events.Define(
        f"n{output_collection[:-1]}", f"return Sum({output_collection}jetmask);"
    )
    for scol in sel_columns:
        # This can be replaced by .Select(Jet_*, ...) when that feature is supported
        events = events.Define(
            output_collection + scol,
            f"return Take({input_collection}{scol}[{output_collection}jetmask], {output_collection}jettake);",
        )
    return events
