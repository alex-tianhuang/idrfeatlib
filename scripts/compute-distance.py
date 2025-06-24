"""
Compute the distance for each iteration of a designs csv.
(see `feature-mimic.py`).
"""


def parse_args():
    import argparse
    parser = argparse.ArgumentParser("compute-distance", description="compute the distance of each iteration in a design csv")
    parser.add_argument("input_sequences", help="whole protein sequences in fasta format")
    parser.add_argument("--input-regions", required=False, help="region boundaries in csv (`ProteinID`, `RegionID`, `Start`, `Stop`) format")
    parser.add_argument("design_csv", help="input designs csv with iterations")
    parser.add_argument("feature_weights_file", help="features csv containing a weight feature vector")
    parser.add_argument("--weights-feature-vector", required=False, default="weights", help="label attached to the weights feature vector")
    parser.add_argument("--feature-file", required=False, help="feature configuration json")
    parser.add_argument("output_csv", help="output csv file")
    return parser.parse_args()

def main():
    args = parse_args()
    from idrfeatlib import FeatureVector
    from idrfeatlib.featurizer import compile_featurizer, Featurizer
    from idrfeatlib.native import compile_native_featurizer
    from idrfeatlib.metric import Metric
    from idrfeatlib.utils import read_fasta, read_nested_csv, read_regions_csv
    import json
    import sys
    import csv
    import math
    import tqdm
    for label, feature_vector in FeatureVector.load(args.feature_weights_file):
        if label == args.weights_feature_vector:
            metric = Metric(feature_vector, feature_vector)
            break
    else:
        raise RuntimeError("could not find feature vector `%s` in %s" % (args.weights_feature_vector, args.feature_weights_file))
    
    if args.feature_file:
        with open(args.feature_file, "r") as file:
            config = json.load(file)
        featurizer, errors = compile_featurizer(config)
    else:
        featurizer, errors = compile_native_featurizer()
    for featname, error in errors.items():
        print("error compiling `%s`: %s" % (featname, error), file=sys.stderr)    
    if featurizer.keys() != metric.weights.as_dict.keys():
        raise RuntimeError("featurizer and metric feature vector `%s` have different features" % args.weights_feature_vector)
    featurizer = Featurizer(featurizer)
    fa = dict(read_fasta(args.input_sequences))
    if args.input_regions is None:
        designs = read_nested_csv(args.design_csv, 2, group_multiple=True)
        with open(args.output_csv, "w") as file:
            w = csv.DictWriter(file, ["ProteinID", "DesignID", "Iteration", "Distance"])
            w.writeheader()
            for protid, ddict in tqdm.tqdm(designs, desc="computing distance"):
                if (entry := fa.get(protid)) is None:
                    continue 
                wt_sequence = entry
                wt_fvec, errors = featurizer.featurize(wt_sequence)
                if len(errors) > 0:
                    for featname, error in errors.items():
                        print("error featurizing wild-type sequence (protid={}) `{}`: {}".format(protid, featname, error), file=sys.stderr)     
                    continue
                for designid, rows in ddict.items():
                    for row in rows:
                        iteration = row["Iteration"]
                        sequence = row["Sequence"]
                        fvec, errors = featurizer.featurize(sequence)
                        if len(errors) > 0:
                            for featname, error in errors.items():
                                print("error featurizing sequence (protid={}, designid={}, iteration={}) `{}`: {}".format(protid, designid, iteration, featname, error), file=sys.stderr)     
                            continue
                        dist = math.sqrt(metric.euclidean_distance_between(wt_fvec, fvec))
                        w.writerow({
                            "ProteinID": protid,
                            "DesignID": designid,
                            "Iteration": iteration,
                            "Distance": dist
                        })
    else:
        regions = read_regions_csv(args.input_regions)
        designs = read_nested_csv(args.design_csv, 3, group_multiple=True)
        with open(args.output_csv, "w") as file:
            w = csv.DictWriter(file, ["ProteinID", "RegionID", "DesignID", "Iteration", "Distance"])
            w.writeheader()
            for protid, rdict in tqdm.tqdm(designs, desc="computing distance"):
                if (entry := fa.get(protid)) is None:
                    continue
                whole_sequence = entry
                assert isinstance(whole_sequence, str)
                if (entryl1 := regions.get(protid)) is None:
                    continue
                for regionid, ddict in rdict.items():
                    if (entryl2 := entryl1.get(regionid)) is None:
                        continue
                    start, stop = entryl2
                    if len((wt_sequence := whole_sequence[start:stop])) != stop - start:
                        print("invalid bounds (protid={},regionid={}) for sequence of length {}: (start={},stop={})".format(protid, regionid, len(whole_sequence), start, stop), file=sys.stderr)
                        continue
                    wt_fvec, errors = featurizer.featurize(wt_sequence)
                    if len(errors) > 0:
                        for featname, error in errors.items():
                            print("error featurizing wild-type sequence (protid={},regionid={}) `{}`: {}".format(protid, regionid,featname, error), file=sys.stderr)     
                        continue
                    for designid, rows in ddict.items():
                        for row in rows:
                            iteration = row["Iteration"]
                            sequence = row["Sequence"]
                            fvec, errors = featurizer.featurize(sequence)
                            if len(errors) > 0:
                                for featname, error in errors.items():
                                    print("error featurizing sequence (protid={}, regionid={}, designid={}, iteration={}) `{}`: {}".format(protid, regionid, designid, iteration, featname, error), file=sys.stderr)     
                                continue
                            dist = math.sqrt(metric.euclidean_distance_between(wt_fvec, fvec))
                            w.writerow({
                                "ProteinID": protid,
                                "RegionID": regionid,
                                "DesignID": designid,
                                "Iteration": iteration,
                                "Distance": dist
                            })
    

    
     
if __name__ == "__main__":
    main()