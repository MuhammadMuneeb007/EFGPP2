"""GA4GH Phenopacket (v2) export of registered phenotypes.

Binary phenotypes with an `ontology_term` become PhenotypicFeatures (excluded = control);
continuous phenotypes become Measurements. When the `phenopackets` package is installed,
every document is round-tripped through its protobuf schema as validation.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from efgpp import __version__
from efgpp.data.registry import Registry, utcnow
from efgpp.project import Project


def _term(definition: dict[str, Any]) -> dict[str, str]:
    term = definition.get("ontology_term") or {}
    return {"id": term.get("id") or f"EFGPP:{definition['phenotype_id']}",
            "label": term.get("label") or definition["name"]}


def build_phenopackets(project: Project) -> dict[str, dict[str, Any]]:
    with Registry.open(project) as reg:
        defs = {d["phenotype_id"]: d for d in reg.rows("SELECT * FROM phenotype_definitions")}
        obs = reg.rows("SELECT * FROM phenotype_observations WHERE NOT is_missing ORDER BY participant_id")
    packets: dict[str, dict[str, Any]] = {}
    for o in obs:
        d = defs[o["phenotype_id"]]
        definition = {**(d.get("definition") or {}), "phenotype_id": d["phenotype_id"], "name": d["name"]}
        pkt = packets.setdefault(o["participant_id"], {
            "id": f"{project.config.project.name}.{o['participant_id']}",
            "subject": {"id": o["participant_id"]},
            "phenotypicFeatures": [], "measurements": [],
            "metaData": {"created": utcnow().isoformat(timespec="seconds") + "Z", "createdBy": f"efgpp {__version__}",
                         "resources": [], "phenopacketSchemaVersion": "2.0"},
        })
        if d["type"] == "binary":
            pkt["phenotypicFeatures"].append({"type": _term(definition), "excluded": o["value_numeric"] == 0.0})
        elif d["type"] == "continuous":
            unit = d.get("units") or "unknown"
            pkt["measurements"].append({"assay": _term(definition), "value": {"quantity": {
                "unit": {"id": f"UCUM:{unit}", "label": unit}, "value": o["value_numeric"]}}})
        else:
            pkt["phenotypicFeatures"].append({"type": _term(definition),
                                              "modifiers": [{"id": f"EFGPP:{o['value_text']}", "label": str(o["value_text"])}]})
    for pkt in packets.values():
        pkt["phenotypicFeatures"] = pkt["phenotypicFeatures"] or None
        pkt["measurements"] = pkt["measurements"] or None
        for k in [k for k, v in pkt.items() if v is None]:
            del pkt[k]
    return packets


def validate_packet(packet: dict[str, Any]) -> str | None:
    try:
        from google.protobuf.json_format import Parse
        from phenopackets import Phenopacket
    except ImportError:
        return None
    try:
        Parse(json.dumps(packet), Phenopacket())
    except Exception as exc:  # noqa: BLE001
        return str(exc)
    return None


def export(project: Project, out_dir: Path | None = None) -> tuple[Path, list[str]]:
    out_dir = out_dir or project.path("reports", "data", "phenopackets")
    out_dir.mkdir(parents=True, exist_ok=True)
    problems = []
    for pid, pkt in build_phenopackets(project).items():
        err = validate_packet(pkt)
        if err:
            problems.append(f"{pid}: {err}")
        (out_dir / f"{pid}.json").write_text(json.dumps(pkt, indent=2), encoding="utf-8")
    return out_dir, problems
