from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl

from efgpp.data.annotation import info_field, lookup_table, write_sites_vcf
from efgpp.data.annotation.opencravat import parse_cravat_tsv
from efgpp.data.annotation.vep import consequence_counts, parse_vep_tab
from efgpp.data.predicted.metaxcan import parse_prediction
from efgpp.data.references.gwas_catalog import sumstats_directory

VARIANTS = pl.DataFrame({"chromosome": ["1", "1", "2"], "position": [100, 200, 300], "variant_id": ["a", "b", "c"],
                         "reference": ["A", "C", "G"], "alternate": ["G", "T", "A"]})


def test_vep_tab_parsing(tmp_path: Path) -> None:
    out = tmp_path / "vep.tsv"
    out.write_text(
        "## VEP output\n"
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tConsequence\tExtra\n"
        "a\t1:100\tG\tENSG1\tENST1\tmissense_variant\tSYMBOL=GENE1;PICK=1\n"
        "a\t1:100\tG\tENSG1\tENST2\tintron_variant\tSYMBOL=GENE1\n"
        "b\t1:200\tT\t-\t-\tintergenic_variant,regulatory_region_variant\tPICK=1\n")
    df = parse_vep_tab(out)
    assert df.height == 3 and "SYMBOL" in df.columns and df.get_column("Gene").to_list()[2] is None
    assert consequence_counts(df) == {"missense_variant": 1, "intergenic_variant": 1, "regulatory_region_variant": 1}


def test_sites_vcf(tmp_path: Path) -> None:
    text = write_sites_vcf(VARIANTS, tmp_path / "s.vcf", "GRCh38").read_text()
    assert "##reference=GRCh38" in text and "1\t100\ta\tA\tG" in text


def test_clinvar_style_lookup(tmp_path: Path) -> None:
    vcf = tmp_path / "clinvar.vcf.gz"
    with gzip.open(vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.1\n##fileDate=20260901\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        fh.write("1\t100\t111\tA\tG,C\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided\n")
        fh.write("chr2\t300\t222\tG\tA\t.\t.\tCLNSIG=Benign\n")
        fh.write("1\t200\t333\tC\tG\t.\t.\tCLNSIG=Benign\n")  # different ALT: must not match
    out = lookup_table(vcf, VARIANTS, select_sql=", ".join(["r.ID AS clinvar_id", info_field("CLNSIG")]))
    got = {r["position"]: r["clnsig"] for r in out.to_dicts()}
    assert got == {100: "Pathogenic", 300: "Benign"}


def test_alphamissense_style_lookup(tmp_path: Path) -> None:
    tsv = tmp_path / "am.tsv.gz"
    with gzip.open(tsv, "wt") as fh:
        fh.write("# Copyright notice\n# more\n#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n")
        fh.write("chr1\t100\tA\tG\thg38\tP1\tT1\tA10G\t0.91\tlikely_pathogenic\n")
    out = lookup_table(tsv, VARIANTS, select_sql="r.am_pathogenicity, r.am_class")
    assert out.row(0) == ("1", 100, "A", "G", "0.91", "likely_pathogenic")


def test_cravat_and_prediction_parsers(tmp_path: Path) -> None:
    rep = tmp_path / "x.variant.tsv"
    rep.write_text("#comment\nchrom\tpos\tref_base\talt_base\tclinvar\n1\t100\tA\tG\tPathogenic\n")
    assert parse_cravat_tsv(rep).row(0) == ("1", "100", "A", "G", "Pathogenic")
    pred = tmp_path / "p.txt"
    pred.write_text("FID\tIID\tENSG1\tENSG2\nS1\tS1\t0.5\t-1.2\n")
    assert parse_prediction(pred).get_column("ENSG2").to_list() == [-1.2]


def test_gwas_catalog_ftp_layout() -> None:
    assert sumstats_directory("GCST90000123").endswith("/GCST90000001-GCST90001000/GCST90000123")
    assert sumstats_directory("GCST000123").endswith("/GCST000001-GCST001000/GCST000123")
