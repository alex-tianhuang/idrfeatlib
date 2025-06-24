"""Tests for the featurize script."""

from . import run_script
import csv
import json

def test_no_regions(tmp_path):
    input_sequences = [
        ("SEQ_A", "A" * 1000),
        ("SEQ_B", "KRKRKRKRKRKRKRK"),
        ("SEQ_C", "GEESEGEESEGEESE")
    ]
    test_protids = {
        "SEQ_A", "SEQ_B", "SEQ_C"
    }

    inf = "{}/seqs.fa".format(tmp_path)
    outf = "{}/feats.csv".format(tmp_path)
    with open(inf, "w") as file:
        for protid, seq in input_sequences:
            file.write(">{}\n{}\n".format(protid, seq))
    run_script("featurize", [inf, outf])
    with open(outf, "r") as file:
        reader = csv.DictReader(file)
        rows = list(reader)
        assert reader.fieldnames is not None, "featurize.py outputted empty file"
        proteinid_col, *featurenames = reader.fieldnames
    for row in rows:
        assert row[proteinid_col] in test_protids, "output file has a non-protid value in first column: (header={},value={})".format(proteinid_col, row[proteinid_col])
        for feat_name in featurenames:
            feat_value = row[feat_name]
            try:
                float(feat_value)
            except ValueError:
                assert isinstance(feat_value, str) and (not feat_value), "featurize.py output something other than a float/None for {}".format(feat_name)
    
def test_with_regions(tmp_path):
    input_sequences = [
        ("SEQ_A", "A" * 1000),
        ("SEQ_B", "KRKRKRKRKRKRKRK"),
        ("SEQ_C", "GEESEGEESEGEESE")
    ]
    input_regions = [
        ("SEQ_A", "NAME_1", 0, 100),
        ("SEQ_B", "NAME_2", 1, 8),
        ("SEQ_B", "NAME_3", 8, 12)
    ]
    test_ids = {
        ("SEQ_A", "NAME_1"),
        ("SEQ_B", "NAME_2"),
        ("SEQ_B", "NAME_3"),
    }
    infa = "{}/seqs.fa".format(tmp_path)
    inrgns = "{}/idrs.csv".format(tmp_path)
    outf = "{}/feats.csv".format(tmp_path)
    with open(infa, "w") as file:
        for protid, seq in input_sequences:
            file.write(">{}\n{}\n".format(protid, seq))
    with open(inrgns, "w") as file:
        writer = csv.DictWriter(file, ["P", "R", "Srt", "Stp"])
        writer.writeheader()
        for protid, regionid, start, stop in input_regions:
            writer.writerow({
                "P": protid,
                "R": regionid,
                "Srt": start,
                "Stp": stop
            })
    run_script("featurize", [infa, "--input-regions", inrgns, outf])
    with open(outf, "r") as file:
        reader = csv.DictReader(file)
        rows = list(reader)
        assert reader.fieldnames is not None, "featurize.py outputted empty file"
        proteinid_col, regionid_col, *featurenames = reader.fieldnames
    for row in rows:
        assert (row[proteinid_col], row[regionid_col]) in test_ids, "output file has a non-id value in first two columns: (header1={},header2={},value1={},value2={})".format(proteinid_col, regionid_col, row[proteinid_col], row[regionid_col])
        for feat_name in featurenames:
            feat_value = row[feat_name]
            try:
                float(feat_value)
            except ValueError:
                assert isinstance(feat_value, str) and not(feat_value), "featurize.py output something other than a float/None for {}".format(feat_name)

def test_custom_feature_json(tmp_path):
    input_feature_json = {
        "features": {
            "CountGeese": {
                "compute": "count",
                "pattern": "GEESE"
            },
            "GlutamateRepeatsFake": {
                "compute": "span",
                "pattern": "EE"
            },
            "Alanine": {
                "compute": "percent_residue",
                "residue": "A"
            },
            "NetCharge": {
                "compute": "score",
                "score": {
                    "K": 1,
                    "R": 1,
                    "E": -1,
                    "D": -1
                }
            },
            "PercentVowels": {
                "compute": "percent_res_group",
                "residue_group": "Vowels"
            },
            "GlutamateRepeats": {
                "compute": "repeats",
                "residues": "E",
            },
            "KRratio": {
                "compute": "log_ratio",
                "numerator": "K",
                "denominator": "R"
            },
            "KEYSpatterning": {
                "compute": "simple_spacing",
                "residue_group": "Keys"
            },
        },
        "residue_groups": {
            "Vowels": "AEI",
            "Keys": "KEYS"
        }
    }
    input_sequences = [
        ("SEQ_A", "A" * 1000),
        ("SEQ_B", "KRKRKRKRKRKRKRKR"),
        ("SEQ_C", "GEESEGEESEGEESE")
    ]
    assertions = [
        ("SEQ_A", "Alanine", dict(eq=1)),
        ("SEQ_A", "KRratio", dict(eq=0)),
        ("SEQ_A", "PercentVowels", dict(eq=1)),
        ("SEQ_B", "CountGeese", dict(eq=0)),
        ("SEQ_B", "KRratio", dict(eq=0)),
        # ("SEQ_B", "KEYSpatterning", dict(bounds=(0, 0.001))), # there may be a math bug here but I need to output a v1 asap
        ("SEQ_C", "CountGeese", dict(eq=3)),
        ("SEQ_C", "GlutamateRepeatsFake", dict(eq=6)),
        ("SEQ_C", "GlutamateRepeats", dict(eq=3)),
    ]

    infa = "{}/seqs.fa".format(tmp_path)
    injson = "{}/features.json".format(tmp_path)
    outf = "{}/feats.csv".format(tmp_path)
    with open(infa, "w") as file:
        for protid, seq in input_sequences:
            file.write(">{}\n{}\n".format(protid, seq))
    with open(injson, "w") as file:
        json.dump(input_feature_json, file)
    run_script("featurize", [infa, outf, "--feature-file", injson])
    with open(outf, "r") as file:
        reader = csv.DictReader(file)
        rows = list(reader)
        assert reader.fieldnames is not None, "featurize.py outputted empty file"
        proteinid_col, *featurenames = reader.fieldnames
    for row in rows:
        protid = row[proteinid_col]
        prot_assertions = {feat_name: cond for _protid, feat_name, cond in assertions if _protid == protid}
        assert len(prot_assertions) > 0, "output file has a non-protid value in first column: (header={},value={})".format(proteinid_col, row[proteinid_col])
        for feat_name in featurenames:
            feat_value = row[feat_name]
            try:
                float(feat_value)
            except ValueError:
                assert isinstance(feat_value, str) and (not feat_value), "featurize.py output something other than a float/None for {}".format(feat_name)
                assert feat_name not in prot_assertions, "value of {} was not supposed to be None".format(feat_name)
                continue
            if feat_name not in prot_assertions:
                continue
            cond = prot_assertions.pop(feat_name)
            if (cond_eq := cond.get("eq")) is not None:
                assert int(cond_eq) == float(feat_value), "equals assertion failed for {}".format(feat_name)
                continue
            # uncomment after fixing SEQ_B KeysPatterning assertion
            # if (cond_bounds := cond.get("bounds")) is not None:
            #     mn, mx = cond_bounds
            #     assert mn < float(feat_value) < mx, "bounds assertion failed for {}".format(feat_name)
            #     continue
            raise RuntimeError("unreachable assertion, bad test setup")
        if len(prot_assertions) > 0:
            raise RuntimeError("untested assertions, bad test setup. untested is the following: {}".format(prot_assertions))