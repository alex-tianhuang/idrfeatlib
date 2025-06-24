"""Tests for the feature mimic script."""
from . import run_script
import csv
import pytest
import json

def common_consts():
    input_sequences = [
        ("COX15", "MLFRNIEVGRQAAKLLTRTSSRLAWQSIGASRNISTIRQQIRKTQLYNFKKTVSIRPFSLSSPVFKPHVASESNPIESRLKTSKNVAYWLIGTSGLVFGIVVLGGLTRLTESGLSITEWKPVTGTLPPMNQKEWEEEFIKYKESPEFKLLNSHIDLDEFKFIFFMEWIHRLWGRAIGAVFILPAVYFAVSKKTSGHVNKRLFGLAGLLGLQGFVGWWMVKSGLDQEQLDARKSKPTVSQ"), # 247 more characters but this test is slow
        ("DED1", "MAELSEQVQNLSINDNNENGYVPPHLRGKPRSARNNSSNYNNNNGGYNGGRGGGSFFSNNRRGGYGNGGFFGGNNGGSRSNGRSGGRWIDGKHVPAPRNEKAEIAIFGVPEDPNFQSSGINFDNYDDIPVDASGKDVPEPITEFTSPPLDGLLLENIKLARFTKPTPVQKYSVPIVANGRDLMACAQTGSGKTGGFLFPVLSESFKTGPSPQPESQGSFYQRKAYPTAVIMAPTRELATQ") # 364 more characters but this test is slow
    ]
    input_regions = [
        ("COX15", "NTERM", 0, 45),
        ("DED1", "nteRM", 0, 90)
    ]
    input_feature_json = dict(
        features={
            "Hydrophobic": dict(compute="percent_res_group", residue_group="Hydrophobic"),
            "ChargePatterning": dict(compute="custom_kappa"),
            "NCPR": dict(compute="score", score={"D": -1, "E": -1, "K": 1, "R": 1}, take_average=True),
            "Glycine": dict(compute="percent_residue", residue="G")
        },
        residue_groups={
            "Hydrophobic": "FYW"
        }
    )
    input_weights = {
        "Hydrophobic": 1,
        "ChargePatterning": 3,
        "NCPR": 2,
        "Glycine": 1
    }
    test_protids = {
        "COX15",
        "DED1"
    }
    test_regions = {
        ("COX15", "NTERM"),
        ("DED1", "nteRM")
    }
    return locals()

def test_syntactically_correct_basic(tmp_path):
    consts = common_consts()
    input_sequences = consts["input_sequences"]
    input_regions = consts["input_regions"]
    input_feature_json = consts["input_feature_json"]
    input_weights = consts["input_weights"]
    test_regions = consts["test_regions"]

    infa = "{}/seqs.fa".format(tmp_path)
    inrgns = "{}/idrs.csv".format(tmp_path)
    injson = "{}/features.json".format(tmp_path)
    inwts = "{}/weights.csv".format(tmp_path)
    outf = "{}/designs.csv".format(tmp_path)
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
    with open(injson, "w") as file:
        json.dump(input_feature_json, file)
    with open(inwts, "w") as file:
        print("Name,{}".format(",".join(input_weights.keys())), file=file)
        print("weights,{}".format(",".join(map(str, input_weights.values()))), file=file)
    run_script("feature-mimic", [infa, "--input-regions", inrgns, inwts, outf, "--feature-file", injson])
    with open(outf, "r") as file:
        reader = csv.DictReader(file)
        rows = list(reader)
        assert reader.fieldnames is not None, "design script output empty file"
        proteinid_col, regionid_col, _designid_col, *other_cols = reader.fieldnames
        other_cols.remove("Sequence")
        other_cols.remove("Time")
        featurenames = other_cols
    assert len(rows) == len(input_regions), "expected {} designs, got the following: {}".format(len(input_regions), rows)
    for row in rows:
        assert (row[proteinid_col], row[regionid_col]) in test_regions, "output file has a non-id value in first two columns: (header1={},header2={},value1={},value2={})".format(proteinid_col, regionid_col, row[proteinid_col], row[regionid_col])
        assert len(row["Sequence"]) > 10, "very short sequence generated: {}".format(row["Sequence"])
        for feat_name in featurenames:
            feat_value = row[feat_name]
            try:
                float(feat_value)
            except ValueError:
                assert isinstance(feat_value, str) and not(feat_value), "feature-mimic.py output something other than a float/None for {}".format(feat_name)

def test_syntactically_correct_with_seeds(tmp_path):
    consts = common_consts()
    input_sequences = consts["input_sequences"]
    input_regions = consts["input_regions"]
    input_feature_json = consts["input_feature_json"]
    input_weights = consts["input_weights"]
    test_regions = consts["test_regions"]

    infa = "{}/seqs.fa".format(tmp_path)
    inrgns = "{}/idrs.csv".format(tmp_path)
    inseed = "{}/seed.csv".format(tmp_path)
    injson = "{}/features.json".format(tmp_path)
    inwts = "{}/weights.csv".format(tmp_path)
    outf = "{}/designs.csv".format(tmp_path)
    SEED = "Seed"
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
    with open(inseed, "w") as file:
        writer = csv.DictWriter(file, ["P", "R", "D", SEED])
        writer.writeheader()
        for protid, regionid, _, _ in input_regions:
            writer.writerow({
                "P": protid,
                "R": regionid,
                "D": 0,
                SEED: 1234567890
            })
    with open(injson, "w") as file:
        json.dump(input_feature_json, file)
    with open(inwts, "w") as file:
        print("Name,{}".format(",".join(input_weights.keys())), file=file)
        print("weights,{}".format(",".join(map(str, input_weights.values()))), file=file)
    run_script("feature-mimic", [infa, "--input-regions", inrgns, inwts, outf, "--feature-file", injson, "--seed", inseed])
    with open(outf, "r") as file:
        reader = csv.DictReader(file)
        rows = list(reader)
        assert reader.fieldnames is not None, "design script output empty file"
        proteinid_col, regionid_col, _designid_col, *other_cols = reader.fieldnames
        other_cols.remove("Sequence")
        other_cols.remove("Time")
        featurenames = other_cols
    assert len(rows) == len(input_regions), "expected {} designs, got the following: {}".format(len(input_regions), rows)
    for row in rows:
        assert (row[proteinid_col], row[regionid_col]) in test_regions, "output file has a non-id value in first two columns: (header1={},header2={},value1={},value2={})".format(proteinid_col, regionid_col, row[proteinid_col], row[regionid_col])
        assert len(row["Sequence"]) > 10, "very short sequence generated: {}".format(row["Sequence"])
        for feat_name in featurenames:
            feat_value = row[feat_name]
            try:
                float(feat_value)
            except ValueError:
                assert isinstance(feat_value, str) and not(feat_value), "feature-mimic.py output something other than a float/None for {}".format(feat_name)

def test_syntactically_correct_without_regions(tmp_path):
    consts = common_consts()
    input_sequences = consts["input_sequences"]
    input_feature_json = consts["input_feature_json"]
    input_weights = consts["input_weights"]
    test_protids = consts["test_protids"]
    
    infa = "{}/seqs.fa".format(tmp_path)
    injson = "{}/features.json".format(tmp_path)
    inwts = "{}/weights.csv".format(tmp_path)
    outf = "{}/designs.csv".format(tmp_path)
    with open(infa, "w") as file:
        for protid, seq in input_sequences:
            file.write(">{}\n{}\n".format(protid, seq))
    with open(injson, "w") as file:
        json.dump(input_feature_json, file)
    with open(inwts, "w") as file:
        print("Name,{}".format(",".join(input_weights.keys())), file=file)
        print("weights,{}".format(",".join(map(str, input_weights.values()))), file=file)
    run_script("feature-mimic", [infa, inwts, outf, "--feature-file", injson])
    with open(outf, "r") as file:
        reader = csv.DictReader(file)
        rows = list(reader)
        assert reader.fieldnames is not None, "design script output empty file"
        proteinid_col, _designid_col, *other_cols = reader.fieldnames
        other_cols.remove("Sequence")
        other_cols.remove("Time")
        featurenames = other_cols
    assert len(rows) == len(input_sequences), "expected {} designs, got the following: {}".format(len(input_sequences), rows)
    for row in rows:
        assert row[proteinid_col] in test_protids, "output file has a non-id value in first column: (header={},value={})".format(proteinid_col, row[proteinid_col])
        assert len(row["Sequence"]) > 10, "very short sequence generated: {}".format(row["Sequence"])
        for feat_name in featurenames:
            feat_value = row[feat_name]
            try:
                float(feat_value)
            except ValueError:
                assert isinstance(feat_value, str) and not(feat_value), "feature-mimic.py output something other than a float/None for {}".format(feat_name)

def test_syntactically_correct_without_regions_with_seed(tmp_path):
    consts = common_consts()
    input_sequences = consts["input_sequences"]
    input_regions = consts["input_regions"]
    input_feature_json = consts["input_feature_json"]
    input_weights = consts["input_weights"]
    test_protids = consts["test_protids"]

    infa = "{}/seqs.fa".format(tmp_path)
    inseed = "{}/seed.csv".format(tmp_path)
    injson = "{}/features.json".format(tmp_path)
    inwts = "{}/weights.csv".format(tmp_path)
    outf = "{}/designs.csv".format(tmp_path)
    SEED = "Seed"
    with open(infa, "w") as file:
        for protid, seq in input_sequences:
            file.write(">{}\n{}\n".format(protid, seq))
    with open(inseed, "w") as file:
        writer = csv.DictWriter(file, ["P", "R", "D", SEED])
        writer.writeheader()
        for protid, regionid, _, _ in input_regions:
            writer.writerow({
                "P": protid,
                "R": regionid,
                "D": 0,
                SEED: 1234567890
            })
    with open(injson, "w") as file:
        json.dump(input_feature_json, file)
    with open(inwts, "w") as file:
        print("Name,{}".format(",".join(input_weights.keys())), file=file)
        print("weights,{}".format(",".join(map(str, input_weights.values()))), file=file)
    run_script("feature-mimic", [infa, inwts, outf, "--feature-file", injson, "--seed", inseed])
    with open(outf, "r") as file:
        reader = csv.DictReader(file)
        rows = list(reader)
        assert reader.fieldnames is not None, "design script output empty file"
        proteinid_col, _designid_col, *other_cols = reader.fieldnames
        other_cols.remove("Sequence")
        other_cols.remove("Time")
        featurenames = other_cols
    assert len(rows) == len(input_sequences), "expected {} designs, got the following: {}".format(len(input_sequences), rows)
    for row in rows:
        assert row[proteinid_col] in test_protids, "output file has a non-id value in first column: (header={},value={})".format(proteinid_col, row[proteinid_col])
        assert len(row["Sequence"]) > 10, "very short sequence generated: {}".format(row["Sequence"])
        for feat_name in featurenames:
            feat_value = row[feat_name]
            try:
                float(feat_value)
            except ValueError:
                assert isinstance(feat_value, str) and not(feat_value), "feature-mimic.py output something other than a float/None for {}".format(feat_name)

def test_syntactically_correct_multiprocessing(tmp_path):
    consts = common_consts()
    input_sequences = consts["input_sequences"]
    input_regions = consts["input_regions"]
    input_feature_json = consts["input_feature_json"]
    input_weights = consts["input_weights"]
    test_regions = consts["test_regions"]

    infa = "{}/seqs.fa".format(tmp_path)
    inrgns = "{}/idrs.csv".format(tmp_path)
    injson = "{}/features.json".format(tmp_path)
    inwts = "{}/weights.csv".format(tmp_path)
    outf = "{}/designs.csv".format(tmp_path)
    N = 10
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
    with open(injson, "w") as file:
        json.dump(input_feature_json, file)
    with open(inwts, "w") as file:
        print("Name,{}".format(",".join(input_weights.keys())), file=file)
        print("weights,{}".format(",".join(map(str, input_weights.values()))), file=file)
    run_script("feature-mimic", [infa, "--input-regions", inrgns, inwts, outf, "--feature-file", injson, "--n-random", str(N), "-np", "2"])
    with open(outf, "r") as file:
        reader = csv.DictReader(file)
        rows = list(reader)
        assert reader.fieldnames is not None, "design script output empty file"
        proteinid_col, regionid_col, _designid_col, *other_cols = reader.fieldnames
        other_cols.remove("Sequence")
        other_cols.remove("Time")
        featurenames = other_cols
    assert len(rows) == len(input_regions) * N, "expected {} designs, got the following: {}".format(len(input_regions) * N, rows)
    for row in rows:
        assert (row[proteinid_col], row[regionid_col]) in test_regions, "output file has a non-id value in first two columns: (header1={},header2={},value1={},value2={})".format(proteinid_col, regionid_col, row[proteinid_col], row[regionid_col])
        assert len(row["Sequence"]) > 10, "very short sequence generated: {}".format(row["Sequence"])
        for feat_name in featurenames:
            feat_value = row[feat_name]
            try:
                float(feat_value)
            except ValueError:
                assert isinstance(feat_value, str) and not(feat_value), "feature-mimic.py output something other than a float/None for {}".format(feat_name)
