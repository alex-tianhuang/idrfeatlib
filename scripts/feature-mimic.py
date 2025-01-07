"""WIP: Script for designing feature-based mimics."""

def parse_args():
    import argparse
    parser = argparse.ArgumentParser("feature-mimic", description="design feature mimics of the given IDR regions")
    inputs = parser.add_argument_group("design-inputs")
    inputs.add_argument("input_sequences", help="whole protein sequences in fasta format")
    inputs.add_argument("feature_weights_file", help="features csv containing a weight feature vector")
    inputs.add_argument("--weights-feature-vector", required=False, default="Label", help="label attached to the weights feature vector")
    inputs.add_argument("--input-regions", required=False, help="region boundaries in csv (`ProteinID`, `RegionID`, `Start`, `Stop`) format")
    parser.add_argument("output_file", help="output csv file with columns (`ProteinID`, `RegionID`, `DesignID`, `Sequence`, ...)")
    rng_seed = parser.add_mutually_exclusive_group()
    rng_seed.add_argument("--n-random", type=int, help="sample this many random query sequences per region")
    rng_seed.add_argument("--seeds-file", help="input csv with (`ProteinId`, `RegionID`, `Seed` | `DesignID`) format")
    parser.add_argument("--feature-file", required=False, help="feature configuration json")
    parser.add_argument("--keep-trajectory", action="store_true", help="when set, save every iteration of the design loop")
    parser.add_argument("--save-seed", action="store_true", help="when set, store the seed in a `Seed` column")
    parser.add_argument("--design-id", required=False, default="{seed}", help="string to format the design id. default is to use the seed")
    return parser.parse_args()

def design_task(designer, query, target, protid, regionid, designid, seed, acceptable_errors=(ArithmeticError, ValueError, KeyError)):
    import sys
    try:
        designer.metric.origin, _ = designer.featurizer.featurize(target, acceptable_errors=())
    except acceptable_errors as e:
        print("cannot featurize target (protid=%s,regionid=%s): %s" % (protid, regionid, e), file=sys.stderr)
        return
    designer.rng.seed(seed)
    AMINOACIDS = list("ACDEFGHIKLMNPQRSTVWY")
    while query is None:
        try_query = "".join(designer.rng.choice(AMINOACIDS) for _ in range(len(target)))
        try:
            designer.featurizer.featurize(try_query, acceptable_errors=())
        except acceptable_errors:
            continue
        query = try_query
    for progress in designer.design_loop(query, acceptable_errors=acceptable_errors):
        ...

def main():
    args = parse_args()
    from benchstuff import Fasta, Regions, RegionsDict, ProteinDict
    from idrfeatlib import FeatureVector
    from idrfeatlib.featurizer import Featurizer, compile_featurizer
    from idrfeatlib.native import compile_native_featurizer
    from idrfeatlib.metric import Metric
    from idrfeatlib.designer import FeatureDesigner
    import json
    import tqdm
    import sys
    import random
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
    if featurizer.keys() != metric.weights.as_dict.keys():
        raise RuntimeError("featurizer and metric feature vector `%s` have different features" % args.weights_feature_vector)
    for featname, error in errors.items():
        print("error compiling `%s`: %s" % (featname, error), file=sys.stderr)    
    designer = FeatureDesigner(featurizer, metric, covergence_threshold=..., good_moves_threshold=..., rng=random.Random())
    LENGTH_THRESHOLD = 30
    SEED_COLNAME = "Seed"
    MAX_SEED = 2 ** 64
    fa = Fasta.load(args.input_sequences)
    fa.assume_unique = True
    tasks = []
    if args.input_regions is None:
        fa = fa.filter(lambda _, seq: len(seq) >= LENGTH_THRESHOLD)
        if args.seeds_file is None:
            n_random = args.n_random or 1
            rng = random.Random()
            seeds = fa.to_protein_dict().map_values(lambda _: [rng.randint(0, MAX_SEED) for _ in range(n_random)])
        else:
            seeds = ProteinDict.load(args.seeds_file, assume_unique=False)
            seeds = seeds.filter(lambda protid, _: protid in fa).map_values(lambda row: row[SEED_COLNAME])
        for protid, prot_seeds in seeds:
            target = fa[protid]
            for counter, seed in enumerate(prot_seeds):
                design_id = args.design_id.format(counter=counter, seed=seed, proteinid=protid)
                tasks.append(
                    (designer, None, target, protid, None, design_id, seed)
                )
    else:
        regions, _ = Regions.load(args.input_regions)
        regions = regions.filter(lambda _p, _r, region: region.len() >= LENGTH_THRESHOLD)
        if args.seeds_file is None:
            n_random = args.n_random or 1
            rng = random.Random()
            seeds = regions.map_values(lambda _: [rng.randint(0, MAX_SEED) for _ in range(n_random)])
        else:
            seeds = RegionsDict.load(args.seeds_file, assume_unique=False)
            seeds = seeds.filter(lambda protid, regionid, _: regionid in (regions.get(protid) or {}))
        for protid, regionid, region_seeds in seeds.iter_nested():
            target_whole = fa[protid]
            start, stop = regions[protid][regionid] 
            target = target_whole[start:stop] # type: ignore
            for counter, seed in enumerate(region_seeds):
                design_id = args.design_id.format(counter=counter, seed=seed, proteinid=protid, regionid=regionid, start=start, stop=stop)
                tasks.append(
                    (designer, None, target, protid, regionid, design_id, seed)
                )
        
raise NotImplementedError("todo finish module")