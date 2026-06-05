"""Config and sample-sheet loader with validation."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd
import yaml


@dataclass
class BarcodeConfig:
    barcode_length: int
    umi_length: int
    ssi_length: int
    ssi_match_length: int

    @property
    def bu_length(self) -> int:
        return self.barcode_length + self.umi_length

    @property
    def total_read_length(self) -> int:
        return self.barcode_length + self.umi_length + self.ssi_length

    @property
    def r2_extract_length(self) -> int:
        return self.umi_length + self.ssi_length


@dataclass
class LibraryFilter:
    enabled: bool
    regex: Optional[str]

    def compile(self) -> Optional[re.Pattern]:
        if not self.enabled or not self.regex:
            return None
        return re.compile(self.regex)


@dataclass
class Bowtie2Config:
    N: int
    L: int
    D: int
    R: int
    i: str
    score_min_local: str
    score_min_global: str
    threads: int
    extra_flags: list[str] = field(default_factory=list)

    def args(self, score_min: str) -> list[str]:
        return [
            "-D", str(self.D), "-R", str(self.R),
            "-N", str(self.N), "-L", str(self.L),
            "-i", self.i,
            "--score-min", score_min,
            "-p", str(self.threads),
            *self.extra_flags,
        ]


@dataclass
class FilteringConfig:
    source_threshold_umi: int
    proj_threshold_umi: int
    cell_body_threshold: Optional[int]


@dataclass
class QualityFilterConfig:
    homopolymer_min_run: int  # 0 or 1 disables; Hyopil default = 7


@dataclass
class NormalizationConfig:
    method: str
    log_factor: float
    spike_method: str = "unique_count"   # "unique_count" (Hyopil normBCmat2) or "umi_sum"


@dataclass
class ClusteringConfig:
    method: str
    metric: str
    num_clusters: int
    cluster_labels: dict[int, str]


@dataclass
class InputConfig:
    raw_fastq_dir: Path
    r1_pattern: str
    r2_pattern: str
    sample_sheet: Path


@dataclass
class ProjectConfig:
    name: str
    output_dir: Path


@dataclass
class Config:
    project: ProjectConfig
    input: InputConfig
    barcode: BarcodeConfig
    library_filter: LibraryFilter
    spike_in_tag: str
    bowtie2: Bowtie2Config
    filtering: FilteringConfig
    quality_filter: QualityFilterConfig
    normalization: NormalizationConfig
    clustering: ClusteringConfig

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Config":
        path = Path(path)
        with path.open("r", encoding="utf-8") as f:
            raw = yaml.safe_load(f)
        base_dir = path.parent
        return cls._from_dict(raw, base_dir)

    @classmethod
    def _from_dict(cls, d: dict, base_dir: Path) -> "Config":
        def resolve(p: str) -> Path:
            pp = Path(p)
            return pp if pp.is_absolute() else (base_dir / pp).resolve()

        return cls(
            project=ProjectConfig(
                name=d["project"]["name"],
                output_dir=resolve(d["project"]["output_dir"]),
            ),
            input=InputConfig(
                raw_fastq_dir=resolve(d["input"]["raw_fastq_dir"]),
                r1_pattern=d["input"]["r1_pattern"],
                r2_pattern=d["input"]["r2_pattern"],
                sample_sheet=resolve(d["input"]["sample_sheet"]),
            ),
            barcode=BarcodeConfig(**d["barcode"]),
            library_filter=LibraryFilter(
                enabled=d["library_filter"].get("enabled", False),
                regex=d["library_filter"].get("regex"),
            ),
            spike_in_tag=d["spike_in"]["tag_sequence"],
            bowtie2=Bowtie2Config(**d["bowtie2"]),
            filtering=FilteringConfig(**d["filtering"]),
            quality_filter=QualityFilterConfig(
                homopolymer_min_run=int(
                    d.get("quality_filter", {}).get("homopolymer_min_run", 7)
                ),
            ),
            normalization=NormalizationConfig(
                method=d["normalization"]["method"],
                log_factor=float(d["normalization"]["log_factor"]),
                spike_method=str(d["normalization"].get("spike_method", "unique_count")),
            ),
            clustering=ClusteringConfig(
                method=d["clustering"]["method"],
                metric=d["clustering"]["metric"],
                num_clusters=d["clustering"]["num_clusters"],
                cluster_labels={int(k): v for k, v in d["clustering"]["cluster_labels"].items()},
            ),
        )


@dataclass
class Sample:
    sample_id: str
    region: str
    side: str
    role: str
    ssi_full_sequence: str
    umi_threshold: int
    sort_order: Optional[int]

    @property
    def region_label(self) -> str:
        return f"{self.region}-{self.side}" if self.side else self.region

    def ssi_match(self, match_length: int) -> str:
        return self.ssi_full_sequence[:match_length]


def load_sample_sheet(path: str | Path) -> list[Sample]:
    df = pd.read_csv(path, dtype={"sample_id": str, "ssi_full_sequence": str})
    required = {"sample_id", "region", "role", "ssi_full_sequence", "umi_threshold"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"sample sheet missing columns: {missing}")

    samples = []
    for _, row in df.iterrows():
        sort_order = row.get("sort_order")
        if pd.isna(sort_order) or sort_order == "":
            sort_order = None
        else:
            sort_order = int(sort_order)
        samples.append(
            Sample(
                sample_id=str(row["sample_id"]),
                region=str(row["region"]),
                side=str(row.get("side", "") or ""),
                role=str(row["role"]).lower(),
                ssi_full_sequence=str(row["ssi_full_sequence"]).upper(),
                umi_threshold=int(row["umi_threshold"]),
                sort_order=sort_order,
            )
        )

    roles = {s.role for s in samples}
    if "source" not in roles:
        raise ValueError("sample sheet must contain at least one row with role=source")
    return samples


def source_samples(samples: list[Sample]) -> list[Sample]:
    return [s for s in samples if s.role == "source"]


def target_samples(samples: list[Sample]) -> list[Sample]:
    return [s for s in samples if s.role == "target"]


def sorted_targets(samples: list[Sample]) -> list[Sample]:
    targets = target_samples(samples)
    keyed = sorted(
        targets,
        key=lambda s: (s.sort_order if s.sort_order is not None else 1_000_000, s.sample_id),
    )
    return keyed
