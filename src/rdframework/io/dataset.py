from __future__ import annotations

from typing import Protocol


class SimpleDatasetProtocol(Protocol):
    def name(self) -> str:
        raise NotImplementedError()

    def xsec(self) -> float:
        raise NotImplementedError()

    def is_mc(self) -> bool:
        raise NotImplementedError()

    def effective_luminosity(self) -> float:
        raise NotImplementedError()

    def latex_name(self) -> str:
        raise NotImplementedError()

    def files(self) -> list[str]:
        raise NotImplementedError()

    def is_skimmed(self) -> bool:
        raise NotImplementedError()

    def skimming_description(self) -> str | None:
        raise NotImplementedError()

    # friendfiles1 = List[str]
    # friendindex1 = str


class SimpleDataset(SimpleDatasetProtocol):
    def __init__(
        self,
        name: str,
        xsec: float,
        is_mc: bool,
        effective_luminosity: float,
        latex_name: str,
        files: list[str] | str,
        is_skimmed: bool,
        skimming_description: str | None,
    ):
        # The name for the dataset, ideally unique within a given datagroup being analyzed
        self._name = name
        self._xsec = xsec
        self._is_mc = is_mc
        self._effective_luminosity = effective_luminosity
        self._latex_name = latex_name
        if isinstance(files, str):
            self._files = [files]
        elif isinstance(files, list):
            self._files = files
        self._is_skimmed = is_skimmed
        if self._is_skimmed and not isinstance(skimming_description, str):
            raise ValueError("skimmed datasets must have a skimming_description")
        else:
            self._skimming_description = skimming_description

    def name(self) -> str:
        return self._name

    def xsec(self) -> float:
        return self._xsec

    def is_mc(self) -> bool:
        return self._is_mc

    def effective_luminosity(self) -> float:
        return self._effective_luminosity

    def latex_name(self) -> str:
        return self._latex_name

    def files(self) -> list[str]:
        return self._files

    def is_skimmed(self) -> bool:
        return self._is_skimmed

    def skimming_description(self) -> str | None:
        return self._skimming_description
