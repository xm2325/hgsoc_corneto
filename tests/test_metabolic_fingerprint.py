from pathlib import Path

from hgsoc_corneto.metabolic.fingerprint import compare_fingerprint, fingerprint_sbml

MINIMAL_SBML = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2">
  <model id="toy">
    <fbc:listOfGeneProducts>
      <fbc:geneProduct fbc:id="G_g1" fbc:label="g1"/>
      <fbc:geneProduct fbc:id="G_g2" fbc:label="g2"/>
    </fbc:listOfGeneProducts>
    <listOfReactions>
      <reaction id="R_biomass_human">
        <fbc:geneProductAssociation>
          <fbc:geneProductRef fbc:geneProduct="G_g1"/>
        </fbc:geneProductAssociation>
      </reaction>
      <reaction id="R_or">
        <fbc:geneProductAssociation><fbc:or>
          <fbc:geneProductRef fbc:geneProduct="G_g1"/>
          <fbc:geneProductRef fbc:geneProduct="G_g2"/>
        </fbc:or></fbc:geneProductAssociation>
      </reaction>
      <reaction id="R_none"/>
    </listOfReactions>
  </model>
</sbml>
"""


def test_fingerprint_sbml(tmp_path: Path) -> None:
    model = tmp_path / "toy.xml"
    model.write_text(MINIMAL_SBML)

    observed = fingerprint_sbml(model)

    assert observed["reactions"] == 3
    assert observed["genes"] == 2
    assert observed["gpr_single_gene"] == 1
    assert observed["gpr_or"] == 1
    assert observed["gpr_no_gene"] == 1
    assert observed["gpr_partition_complete"] is True
    assert observed["biomass_sbml_ids"] == ["R_biomass_human"]
    assert compare_fingerprint(observed, {"reactions": 3, "genes": 2}) == {}


def test_compare_fingerprint_reports_only_mismatches() -> None:
    assert compare_fingerprint({"reactions": 3}, {"reactions": 4, "genes": 2}) == {
        "reactions": {"expected": 4, "observed": 3},
        "genes": {"expected": 2, "observed": None},
    }
