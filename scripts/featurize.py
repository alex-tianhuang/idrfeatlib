"""
Script for computing the proteome feature vectors.

Example
-------
`$ python featurize.py input.fasta output.csv`
"""

def parse_args():
    import argparse
    parser = argparse.ArgumentParser("featurize", description="compute the feature vectors of a fasta file")
    parser.add_argument("input_sequences", help="input fasta file")
    parser.add_argument("--input-regions", required=False, help="input regions csv, format (`ProteinID`, `RegionID`, `Start`, `Stop`)")
    parser.add_argument("output_file", help="output feature file (will be a csv)")
    parser.add_argument("--feature-file", required=False, help="feature configuration json")
    return parser.parse_args()

def main():
    args = parse_args()
    from benchstuff import Fasta
    from benchstuff import Regions
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
    fa = Fasta.load(args.input_sequences)
    if args.input_regions is None:
        seqs = tqdm.tqdm(fa, total=len(fa), desc="featurizing sequences")
        feature_vectors, errors = featurizer.featurize_to_matrices(seqs)
        for proteinid, error in errors.items():
            print("error for `%s`: %s" % (proteinid, error), file=sys.stderr)
        FeatureVector.dump(list(feature_vectors.items()), args.output_file, "ProteinID")
    else:
        regions, _ = Regions.load(args.input_regions)
        Fasta.assume_unique = True
        seqs = []
        for protid, regionid, (start, stop) in regions.iter_nested():
            if (entry := fa.get(protid)) is None:
                continue
            _, whole_seq = entry
            assert isinstance(whole_seq, str), type(whole_seq)
            seq = whole_seq[start:stop]
            if len(seq) != stop - start:
                print("invalid region `%s` for protein `%s` (start=%d,stop=%d,seqlen=%s)" % (regionid, protid, start, stop, len(whole_seq)))
                continue
            seqs.append(((protid, regionid), seq))
        feature_vectors, errors = featurizer.featurize_to_matrices(tqdm.tqdm(seqs, desc="featurizing sequences"))
        for (protid, regionid), error in errors.items():
            print("error for (protid=%s, regionid=%s): %s" % (protid, regionid, error), file=sys.stderr)
        FeatureVector.dump(list(feature_vectors.items()), args.output_file, ("ProteinID", "RegionID")) 


    
     
if __name__ == "__main__":
    main()