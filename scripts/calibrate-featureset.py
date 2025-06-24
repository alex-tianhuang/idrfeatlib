"""
Script that fills in the expected motif frequencies + the expected aa frequencies,
then puts `subtract_expected` on all "count" and "repeat" features.

Example
-------
`$ python calibrate-featureset.py input-sequences.fasta default output-featureset.json`
"""

def parse_args():
    import argparse
    parser = argparse.ArgumentParser("calibrate-featureset", description="compute the expected motif frequencies + the expected aa frequencies")
    parser.add_argument("input_sequences", help="input fasta file")
    parser.add_argument("--input-regions", required=False, help="input regions csv, format (`ProteinID`, `RegionID`, `Start`, `Stop`)")
    parser.add_argument("input_features", help="input feature config file (will be a json). type in `default` for the default feature file")
    parser.add_argument("output_features", help="output feature config file (will be a json)")
    return parser.parse_args()

def main(args):
    from idrfeatlib.utils import read_regions_csv, read_fasta, iter_nested
    from idrfeatlib.featurizer import Featurizer, compile_featurizer
    from idrfeatlib.native import compile_native_featurizer
    import idrfeatlib
    import json
    import tqdm
    import sys
    import os
    if args.input_features != "default":
        with open(args.input_features, "r") as file:
            config = json.load(file)
        featurizer, errors = compile_featurizer(config)
    else:
        with open(
            os.path.join(os.path.dirname(idrfeatlib.__file__), "native-features.json"), "r"
        ) as file:
            config = json.load(file)
        featurizer, errors = compile_native_featurizer(config)
    for featname, error in errors.items():
        print("error compiling `%s`: %s" % (featname, error), file=sys.stderr)
    for featname, feature in config["features"].items():
        if (compute := feature.get("compute")) is None:
            continue
        if compute != "count":
            _ = featurizer.pop(featname, None)
    featurizer = Featurizer(featurizer)
    fa = dict(read_fasta(args.input_sequences))
    if args.input_regions is None:
        seqs = fa.items()
        feature_vectors, errors = featurizer.featurize_to_matrices(tqdm.tqdm(seqs, total=len(seqs), desc="featurizing sequences"))
        for proteinid, error in errors.items():
            print("error for `%s`: %s" % (proteinid, error), file=sys.stderr)
        
    else:
        regions = read_regions_csv(args.input_regions)
        seqs = []
        for protid, regionid, (start, stop) in iter_nested(regions, 2):
            if (entry := fa.get(protid)) is None:
                continue
            whole_seq = entry
            assert isinstance(whole_seq, str), type(whole_seq)
            seq = whole_seq[start:stop]
            if len(seq) != stop - start:
                print("invalid region `%s` for protein `%s` (start=%d,stop=%d,seqlen=%s)" % (regionid, protid, start, stop, len(whole_seq)))
                continue
            seqs.append(((protid, regionid), seq))
        feature_vectors, errors = featurizer.featurize_to_matrices(tqdm.tqdm(seqs, desc="featurizing sequences"))
        for (protid, regionid), error in errors.items():
            print("error for (protid=%s, regionid=%s): %s" % (protid, regionid, error), file=sys.stderr)
    seqs = dict(seqs)
    counting_features = []
    repeat_features = []
    config["motif_frequencies"] = motif_frequencies = config.get("motif_frequencies") or {}
    for featname, feature in config["features"].items():
        if (compute := feature.get("compute")) is None:
            continue
        if compute == "repeat":
            repeat_features.append(feature)
            continue
        if compute != "count":
            continue
        if feature.get("subtract_expected"):
            continue
        if (pattern := feature.get("pattern")) is None:
            continue
        counting_features.append((featname, feature))
    running_sums = [0.0 for _ in counting_features]
    running_counts = [0 for _ in counting_features]
    AMINOACIDS = {aa: i for i, aa in enumerate("ACDEFGHIKLMNPQRSTVWY")}
    running_composition_counts = [0 for _ in AMINOACIDS]
    running_aa_count = 0
    for seq_id, seq in seqs.items():
        fvec = feature_vectors[seq_id]
        running_aa_count += len(seq)
        for ch in seq:
            if (i := AMINOACIDS.get(ch)) is None:
                running_aa_count -= 1
                continue
            running_composition_counts[i] += 1
        for i, (featname, feature) in enumerate(counting_features):
            featvalue = fvec.as_dict[featname]
            if feature.get("average") or feature.get("take_average"):
                count = featvalue * len(seq)
            else:
                count = featvalue
            pattern = feature["pattern"]
            if (denominator := len(seq) - len(pattern)) <= 0:
                continue
            running_sums[i] += count / denominator
            running_counts[i] += 1
    config["residue_frequencies"] = {
        aa: running_composition_counts[i] / running_aa_count for aa, i in AMINOACIDS.items()
    }
    for i, (featname, feature) in enumerate(counting_features):
        running_sum = running_sums[i]
        running_count = running_counts[i]
        if running_count == 0:
            continue
        pattern = feature["pattern"]
        motif_frequencies[pattern] = running_sum / running_count
        feature["subtract_expected"] = True
    for feature in repeat_features:
        feature["subtract_expected"] = True
    with open(args.output_features, "w") as file:
        json.dump(config,file)
    
     
if __name__ == "__main__":
    main(parse_args())