"""
Script for computing the mean and variance of features over a fasta file.

Example
-------
`$ python compute-metric.py input-sequences.fasta output.csv`
`$ python compute-metric.py input-features.csv output.csv`
"""

def parse_args():
    import argparse
    parser = argparse.ArgumentParser("compute-metric", description="compute the mean and weights of a given feature vector spreadsheet")
    parser.add_argument("input_file", help="input features csv file")
    parser.add_argument("output_file", help="output features csv file")
    parser.add_argument("--primary-id", required=False, help="column name containing the fasta headers in the INPUT file")
    parser.add_argument("--output-label", required=False, default="Label", help="column labelling the origin or weights vector in the OUTPUT file")
    parser.add_argument("--feature-file", required=False, help="feature configuration json")
    return parser.parse_args()

def main():
    args = parse_args()
    from idrfeatlib import FeatureVector
    from idrfeatlib.featurizer import Featurizer, compile_featurizer
    from idrfeatlib.native import compile_native_featurizer
    from idrfeatlib.metric import Metric
    import json
    import sys
    import math
    import tqdm
    if args.feature_file:
        with open(args.feature_file, "r") as file:
            config = json.load(file)
        featurizer, errors = compile_featurizer(config)
    else:
        featurizer, errors = compile_native_featurizer()
    for featname, error in errors.items():
        print("error compiling `%s`: %s" % (featname, error), file=sys.stderr)        
    featurizer = Featurizer(featurizer)
    feature_vectors = (fvec for _, fvec in FeatureVector.load(args.input_file))
    feature_vectors = tqdm.tqdm(feature_vectors, desc="computing mean + var")
    _, mean, variance = FeatureVector.cmv(feature_vectors)
    weights = variance.map_values(lambda x: 1 / math.sqrt(x) if x > 0 else 0)
    Metric(mean, weights).dump(args.output_file, primary_id=args.output_label)

if __name__ == "__main__":
    main()