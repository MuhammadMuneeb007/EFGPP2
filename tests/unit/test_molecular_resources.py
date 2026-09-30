"""Model resource parsing: SpliceAI output, OmicsPred scoring/validation files, Zenodo checks, ancestry."""

from __future__ import annotations

from pathlib import Path

import httpx
import pytest

from efgpp.data.annotation.spliceai import parse_spliceai_info, parse_spliceai_vcf
from efgpp.data.predicted.engine import ancestry_match_status
from efgpp.data.predicted.models import read_scoring_file
from efgpp.data.references.molecular import omicspred_rows, read_validation
from efgpp.resources.downloader import DownloadError
from efgpp.resources.zenodo import fetch_file, record_file

# Header and rows as published by OmicsPred (scoring_files_hm_38)
SCORING = """##OMICSPRED SCORE INFORMATION
#omicspred_id=OPGS003420
#pgs_name=sldlc
#trait_type=metabolomics
#trait_reported=Cholesterol in small LDL
#genome_build=GRCh37
##HARMONIZATION DETAILS
#HmPOS_build=GRCh38
rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_source\thm_rsID\thm_chr\thm_pos\thm_inferOtherAllele
rs12732125\t1\t55470153\tC\tT\t0.0437730082617823\tENSEMBL\trs12732125\t1\t55004480\t
rs3976734\t1\t55489960\tA\tG\t0.0132471697744438\tENSEMBL\trs3976734\t1\t55024287\t
"""


def test_spliceai_scores_kept_raw() -> None:
    rows = parse_spliceai_info("G|GENE1|0.01|0.02|0.91|0.00|-2|10|3|-5,G|GENE2|0.30|.|0.10|0.05|1|2|3|4")
    assert rows[0]["spliceai_max_score"] == 0.91 and rows[0]["spliceai_effect"] == "donor_gain"
    assert rows[0]["DP_DG"] == 3 and rows[1]["DS_AL"] is None and rows[1]["spliceai_effect"] == "acceptor_gain"


def test_spliceai_vcf(tmp_path: Path) -> None:
    vcf = tmp_path / "s.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "1\t100\t1:100:A:G\tA\tG\t.\t.\tSpliceAI=G|G1|0.2|0.0|0.0|0.1|1|2|3|4\n"
                   "1\t200\t1:200:C:T\tC\tT\t.\t.\t.\n", encoding="utf-8")
    t = parse_spliceai_vcf(vcf)
    assert t.height == 1 and t.row(0, named=True)["variant_key"] == "1:100:A:G"


def test_omicspred_scoring_file_uses_harmonized_grch38(tmp_path: Path) -> None:
    f = tmp_path / "OPGS003420_hmPOS_GRCh38.txt"
    f.write_text(SCORING, encoding="utf-8")
    header, w = read_scoring_file(f)
    assert header["_genome_build"] == "GRCh38" and header["omicspred_id"] == "OPGS003420"
    assert w.get_column("position").to_list() == [55004480, 55024287]  # hm_pos, not the GRCh37 position
    assert w.select("effect_allele", "other_allele").row(0) == ("C", "T")


def test_omicspred_validation_and_filters(tmp_path: Path) -> None:
    v = tmp_path / "val"
    v.write_text("OMICSPRED ID\tTrait ID\tInternal_R2\tUKB_R2\nOPGS003420\tsldlc\t0.129\t0.064\n", encoding="utf-8")
    t = read_validation(v)
    assert t is not None and t.row(0) == ("OPGS003420", 0.129, "Internal")
    cat = {"results": [
        {"id": "OPD000003", "name": "INTERVAL Nightingale", "omics_type": "metabolite",
         "platform": {"name": "Nightingale"}, "tissue": {"label": "blood serum"}, "scores_count": 141,
         "samples_training": [{"sample_number": 37359, "ancestry_broad": "European",
                               "cohorts": [{"name_short": "INTERVAL"}]}]},
        {"id": "OPD000001", "name": "INTERVAL SomaScan", "omics_type": "protein", "platform": {"name": "Somalogic"},
         "tissue": {"label": "blood plasma"}, "scores_count": 2384, "samples_training": []}]}
    assert [r["id"] for r in omicspred_rows(cat, modality="metabolomics")] == ["OPD000003"]
    assert [r["id"] for r in omicspred_rows(cat, modality="proteomics", platform="soma")] == ["OPD000001"]
    assert omicspred_rows(cat, ancestry="european")[0]["training_n"] == 37359


def test_ancestry_match_status() -> None:
    assert ancestry_match_status("European", {"EUR"}) == "matched"
    assert ancestry_match_status("European, East Asian", {"EUR", "EAS"}) == "mixed"
    assert ancestry_match_status("European", {"EUR", "AFR"}) == "partial"
    assert ancestry_match_status("European", {"AFR"}) == "mismatched"
    assert ancestry_match_status(None, {"EUR"}) == "unknown" and ancestry_match_status("European", set()) == "unknown"


def _zenodo_client(md5: str, content: bytes) -> httpx.Client:
    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.path == "/api/records/3518299":
            return httpx.Response(200, json={"id": 3518299, "doi": "10.5281/zenodo.3518299",
                                             "metadata": {"title": "t", "license": {"id": "cc-by-4.0"}},
                                             "files": [{"key": "mashr_eqtl.tar", "size": len(content),
                                                        "checksum": f"md5:{md5}",
                                                        "links": {"self": "https://zenodo.test/f"}}]})
        return httpx.Response(200, content=content)

    return httpx.Client(transport=httpx.MockTransport(handler), base_url="https://zenodo.org")


def test_zenodo_checksum_verified_and_changes_refused(tmp_path: Path) -> None:
    import hashlib

    content = b"model archive"
    md5 = hashlib.md5(content).hexdigest()  # noqa: S324
    with _zenodo_client(md5, content) as c:
        assert record_file("3518299", "mashr_eqtl.tar", c).md5 == md5
        path, sha, zf = fetch_file("3518299", "mashr_eqtl.tar", tmp_path, client=c, expected_md5=md5)
        assert path.read_bytes() == content and sha == hashlib.sha256(content).hexdigest()
        assert zf.manifest()["source_checksum"] == f"md5:{md5}"
        with pytest.raises(DownloadError, match="upstream file changed"):
            fetch_file("3518299", "mashr_eqtl.tar", tmp_path / "x", client=c, expected_md5="0" * 32)
    with _zenodo_client("f" * 32, content) as c, pytest.raises(DownloadError, match="MD5 mismatch"):
        fetch_file("3518299", "mashr_eqtl.tar", tmp_path / "y", client=c)
