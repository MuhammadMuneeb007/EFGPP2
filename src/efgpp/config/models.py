"""Shared pydantic building blocks for EFGPP configuration files."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, TypeVar

import yaml
from pydantic import BaseModel, ConfigDict

M = TypeVar("M", bound="StrictModel")


class StrictModel(BaseModel):
    """Configuration model that rejects unknown keys so typos surface immediately."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True, use_enum_values=False)

    def to_yaml_dict(self) -> dict[str, Any]:
        return self.model_dump(mode="json", exclude_none=True)

    def checksum(self) -> str:
        payload = json.dumps(self.model_dump(mode="json"), sort_keys=True, default=str)
        return hashlib.sha256(payload.encode()).hexdigest()


class OpenModel(StrictModel):
    """Model that tolerates provider-specific extra keys (e.g. in resources.yaml)."""

    model_config = ConfigDict(extra="allow", validate_assignment=True)


def read_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    return data or {}


def write_yaml(path: Path, data: dict[str, Any], header: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = yaml.safe_dump(data, sort_keys=False, allow_unicode=True, default_flow_style=False)
    if header:
        text = "".join(f"# {line}\n" if line else "#\n" for line in header.splitlines()) + "\n" + text
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def load_model(cls: type[M], path: Path) -> M:
    return cls.model_validate(read_yaml(path))


def file_checksum(path: Path) -> str | None:
    if not path.exists():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()
