"""
Script for computing the proteome feature vectors.

Example
-------
`$ python featurize.py input.fasta output.csv`
"""

def parse_args():
    import argparse
    parser = argparse.ArgumentParser("featurize", description="compute the feature vectors of a fasta file")
    parser.add_argument("input_file", help="input fasta file")
    parser.add_argument("output_file", help="output feature file (will be a csv)")
    parser.add_argument("--primary-id", required=False, default="ProteinID", help="column name containing the fasta headers")
    parser.add_argument("--feature-file", required=False, help="feature configuration json")
    return parser.parse_args()

def main():
    args = parse_args()
    from benchstuff import Fasta
    from idrfeatlib import FeatureVector
    from idrfeatlib.featurizer import Featurizer, compile_featurizer
    from idrfeatlib.native import compile_native_featurizer
    import json
    import tqdm
    import sys
    if args.feature_file:
        with open(args.feature_file, "r") as file:
            config = json.load(file)
        featurizer, errors = compile_featurizer(config)
    else:
        featurizer, errors = compile_native_featurizer()
    for featname, error in errors.items():
        print("error compiling `%s`: %s" % (featname, error), file=sys.stderr)        
    featurizer = Featurizer(featurizer)
    fa = Fasta.load(args.input_file)
    seqs = tqdm.tqdm(fa, total=len(fa), desc="featurizing sequences")
    feature_vectors, errors = featurizer.featurize_to_matrices(seqs)
    for proteinid, error in errors.items():
        print("error for `%s`: %s" % (proteinid, error), file=sys.stderr)
    FeatureVector.dump(list(feature_vectors.items()), args.output_file, args.primary_id)
     
if __name__ == "__main__":
    main()