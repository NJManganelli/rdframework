from __future__ import annotations

from typing import Any


def vidUnpackedWP(
    events: Any, return_columns: bool = True, input_collection: str = "Electron_"
) -> Any:
    """Return dataframe with columns of the cuts in the electron cutBasedID,
    e.g. Electron_GsfEleEInverseMinusPInverseCut will be 0 (fail), 1, 2, 3, or 4 (tight)"""
    columns = events.GetColumnNames()
    ret_cols = []
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
        # hanky, but more consistent with coffea implementation line-by-line
        ret_cols.append(name)
        # Replace with Redefine when it can be used to 'define if undefined, overwrite if defined'
        if f"{input_collection}{name}" in columns:
            continue
        events = events.Define(
            f"{input_collection}{name}",
            (
                f"ROOT::VecOps::RVec<int> tmp; for(auto v : {input_collection}vidNestedWPBitmap){{\n"
                f"  tmp.push_back((v >> {shift}) & 0b111);\n"
                "} //end loop over vidNestedWPBitmap\n"
                "return tmp;"
            ),
        )
    return events, ret_cols


def vidUnpackedWPSelection(electrons: Any, level: int) -> Any:
    """Return a dictionary of boolean masks for the electron cutBasedID,
    e.g. results["GsfEleEInverseMinusPInverseCut"] will be True if the result value is >= level"""
    raise NotImplementedError()
    # results = dict()
    # for name, cut_level in vidUnpackedWP(electrons).items():
    #     results[name] = cut_level >= level
    # return results


def select_electrons_cutBased(
    events: Any,
    input_collection: str,
    output_collection: str,
    columns: list[str] | None,
    electron_id: int | str,
    min_pt: float,
    max_eta: float = 2.5,
    max_ip3d_barrel: float = 0.05,
    max_ip3d_endcap: float = 0.10,
    max_dz_barrel: float = 0.1,
    max_dz_endcap: float = 0.2,
    invert_cuts: list[str] | None = None,
    sort_column: str | None = "pt",
    sort_reverse: bool = True,
) -> Any:
    if not input_collection.endswith("_"):
        input_collection += "_"
    if not output_collection.endswith("_"):
        output_collection += "_"
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
    mask = f"""auto emask =
    (
        (
            (abs({input_collection}eta) < 1.4442)
            && (abs({input_collection}ip3d) < {max_ip3d_barrel})
            && (abs({input_collection}dz) < {max_dz_barrel})
        )
        ||
        (
            (abs({input_collection}eta) > 1.5660)
            && (abs({input_collection}eta) <= {max_eta})
            && (abs({input_collection}ip3d) < {max_ip3d_endcap})
            && (abs({input_collection}dz) < {max_dz_endcap})
        )
    )
    &&
    ({input_collection}pt >= {min_pt})"""
    if isinstance(invert_cuts, list) and len(invert_cuts) > 0:
        cln_invert_cuts = [cut.split("_")[-1] for cut in invert_cuts]
        # add unpacked columns to dataset
        events, vid_cuts = vidUnpackedWP(
            events, return_columns=True, input_collection=input_collection
        )
        for name in vid_cuts:
            if name not in cln_invert_cuts:
                # Need the cut to have passed at the minimum level
                mask += f"\n && ({input_collection}{name} >= {e_id})"
            else:
                # Unfortunately, not all are evaluated for all levels, which makes things confusing
                mask += f"\n && !({input_collection}{name} >= {e_id})"
    else:
        mask += f"\n && ({input_collection}cutBased >= {e_id})"
    mask += "; return emask;"

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
    events = events.Define(f"{output_collection}elmask", mask)

    # Now define all our columns... what about sorting!?
    if isinstance(columns, list):
        sel_columns = [
            col[len(input_collection) :]
            for col in columns
            if col.startswith(input_collection)
        ]
    elif isinstance(columns, str):
        raise NotImplementedError("regexp not currently supported")
    else:
        sel_columns = [col[len(input_collection) :] for col in avail_columns] + [
            f"{input_collection}idx"
        ]

    if sort_column:
        # Build a take vector that must be applied to the SLICED variables
        # e.g. Take(Electron_eta[ele_mask], sort_indices_for_slice)
        if sort_reverse:
            events = events.Define(
                f"{output_collection}eltake",
                f"return Reverse(Argsort({input_collection}{sort_column}[{output_collection}elmask]));",
            )
        else:
            events = events.Define(
                f"{output_collection}eltake",
                f"return Argsort({input_collection}{sort_column}[{output_collection}elmask]);",
            )
    else:
        events = events.Define(
            f"{output_collection}eltake",
            f"return {input_collection}idx[{output_collection}elmask];",
        )

    events = events.Define(
        f"n{output_collection[:-1]}", f"return Sum({output_collection}elmask);"
    )
    for scol in sel_columns:
        # This can be replaced by .Select(Jet_*, ...) when that feature is supported
        events = events.Define(
            output_collection + scol,
            f"return Take({input_collection}{scol}[{output_collection}elmask], {output_collection}eltake);",
        )
    return events


def select_muons_cutBased(
    events: Any,
    input_collection: str,
    output_collection: str,
    columns: list[str] | None,
    muon_id: str,
    muon_iso: int | str,
    min_pt: float,
    max_eta: float = 2.4,
    max_ip3d: float = 0.10,
    max_dz: float = 0.2,
    invert_iso: bool = False,
    sort_column: str | None = "pt",
    sort_reverse: bool = True,
) -> Any:
    if not input_collection.endswith("_"):
        input_collection += "_"
    if not output_collection.endswith("_"):
        output_collection += "_"

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
    mask = f"""auto mmask = (
        ({input_collection}pt >= {min_pt})
        && (abs({input_collection}eta) <= {max_eta})
        && (abs({input_collection}ip3d) < {max_ip3d})
        && (abs({input_collection}dz) < {max_dz})
    )"""
    if invert_iso:
        mask += f"\n && ({input_collection}pfIsoId < {m_iso})"
    else:
        mask += f"\n && ({input_collection}pfIsoId >= {m_iso})"

    # apply the cutbased ID, which is stored in separate boolean branches/fields
    if muon_id.lower() == "fail":
        pass  # accept all the muons...
    elif muon_id.lower() == "veto":
        raise ValueError("Muons do not have a veto cutBased working point")
    elif muon_id.lower() == "loose":
        mask += f"\n && ({input_collection}looseId)"
    elif muon_id.lower() == "medium":
        mask += f"\n && ({input_collection}mediumId)"
    elif muon_id.lower() == "tight":
        mask += f"\n && ({input_collection}tightId)"

    mask += "; return mmask;"

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
    events = events.Define(f"{output_collection}mumask", mask)

    # Now define all our columns... what about sorting!?
    if isinstance(columns, list):
        sel_columns = [
            col[len(input_collection) :]
            for col in columns
            if col.startswith(input_collection)
        ]
    elif isinstance(columns, str):
        raise NotImplementedError("regexp not currently supported")
    else:
        sel_columns = [col[len(input_collection) :] for col in avail_columns] + [
            f"{input_collection}idx"
        ]

    if sort_column:
        # Build a take vector that must be applied to the SLICED variables
        # e.g. Take(Muon_eta[mu_mask], sort_indices_for_slice)
        if sort_reverse:
            events = events.Define(
                f"{output_collection}mutake",
                f"return Reverse(Argsort({input_collection}{sort_column}[{output_collection}mumask]));",
            )
        else:
            events = events.Define(
                f"{output_collection}mutake",
                f"return Argsort({input_collection}{sort_column}[{output_collection}mumask]);",
            )
    else:
        events = events.Define(
            f"{output_collection}mutake",
            f"return {input_collection}idx[{output_collection}mumask];",
        )
    events = events.Define(
        f"n{output_collection[:-1]}", f"return Sum({output_collection}mumask);"
    )
    for scol in sel_columns:
        # This can be replaced by .Select(Jet_*, ...) when that feature is supported
        events = events.Define(
            output_collection + scol,
            f"return Take({input_collection}{scol}[{output_collection}mumask], {output_collection}mutake);",
        )
    return events
