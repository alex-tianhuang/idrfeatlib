"""
Script for finding sequences with similar features as those in wild-type sequences.

Example
-------
`$ python feature-mimic.py input.fasta --input-regions idr-bounds.csv output.csv`
"""

def parse_args():
    import argparse
    parser = argparse.ArgumentParser("feature-mimic", description="design feature mimics of the given IDR regions")
    inputs = parser.add_argument_group("design-inputs")
    inputs.add_argument("input_sequences", help="whole protein sequences in fasta format")
    inputs.add_argument("feature_weights_file", help="features csv containing a weight feature vector")
    inputs.add_argument("--weights-feature-vector", required=False, default="weights", help="label attached to the weights feature vector")
    inputs.add_argument("--input-regions", required=False, help="region boundaries in csv (`ProteinID`, `RegionID`, `Start`, `Stop`) format")
    parser.add_argument("output_file", help="output csv file with columns (`ProteinID`, `RegionID`, `DesignID`, `Sequence`, ...)")
    seed_args = parser.add_mutually_exclusive_group()
    seed_args.add_argument("--n-random", type=int, help="sample this many random query sequences per region")
    seed_args.add_argument("--seeds-file", help="input csv with (`ProteinID`, `RegionID`, `Seed` | `DesignID`) format")
    seed_args.add_argument("--query-file", help="csv full of query sequences (in column `Sequence`), labelled by `ProteinID`, `DesignID` (and `RegionID` if applicable)")
    parser.add_argument("--feature-file", required=False, help="feature configuration json")
    parser.add_argument("--keep-trajectory", action="store_true", help="when set, save every iteration of the design loop")
    parser.add_argument("--save-seed", action="store_true", help="when set, store the seed in a `Seed` column")
    parser.add_argument("--design-id", required=False, default="{counter}", help="string to format the design id. default is to use a counter")
    parser.add_argument("-np", "--n-processes", type=int, required=False, default=1, help="number of processes. requires libraries: pathos + tqdm_pathos")
    parser.add_argument("--greedy", action="store_true", help="use a greedy optimization algorithm instead of the default fast one.")
    return parser.parse_args()

def init_subprocess(stderr_lock, output_lock):
    global STDERR_LOCK
    STDERR_LOCK = stderr_lock
    global OUTPUT_LOCK
    OUTPUT_LOCK = output_lock

def design_task(query, target, protid, regionid, designid, seed, designer, colnames, args, acceptable_errors=(ArithmeticError, ValueError, KeyError)):
    import sys
    import csv
    from idrfeatlib.featurizer import Featurizer
    featurizer = Featurizer(designer.featurizer)
    try:
        designer.metric.origin, _ = featurizer.featurize(target, acceptable_errors=())
    except acceptable_errors as e:
        with STDERR_LOCK:
            print("cannot featurize target (protid=%s,regionid=%s): %s" % (protid, regionid, e), file=sys.stderr)
        return
    designer.rng.seed(seed)
    AMINOACIDS = list("ACDEFGHIKLMNPQRSTVWY")
    MAX_RETRIES = 15
    if query is None:
        for _ in range(MAX_RETRIES):
            try_query = "".join(designer.rng.choice(AMINOACIDS) for _ in range(len(target)))
            try:
                featurizer.featurize(try_query, acceptable_errors=())
            except acceptable_errors:
                continue
            query = try_query
            break
        else:
            with STDERR_LOCK:
                print("cannot generate query with all features (protid=%s,regionid=%s,length=%d,seed=%d)" % (protid, regionid, len(target), seed), file=sys.stderr)
            return
    try:
        if args.keep_trajectory:
            save = []
            for progress in designer.design_loop(query, acceptable_errors=acceptable_errors):
                save.append(progress)
        else:
            for progress in designer.design_loop(query, acceptable_errors=acceptable_errors):
                progress.pop("Iteration") 
            save = [progress]
    except acceptable_errors:
        with STDERR_LOCK:
            print("query did not have all features (protid=%s,regionid=%s,seed=%d)" % (protid, regionid, seed), file=sys.stderr)
        return
    with OUTPUT_LOCK:
        with open(args.output_file, "a") as file:
            writer = csv.DictWriter(file, colnames)
            for row in save:
                row["ProteinID"] = protid
                row["DesignID"] = designid
                if regionid is not None:
                    row["RegionID"] = regionid
                if args.save_seed:
                    row["Seed"] = seed
                writer.writerow(row)
def design_all(num_processes, tasks):
    if num_processes > 1:
        import multiprocessing
        import pathos # type: ignore
        import tqdm_pathos # type: ignore
        stderr_lock = multiprocessing.Lock()
        output_lock = multiprocessing.Lock()
        pool = pathos.multiprocessing.Pool(num_processes, initializer=init_subprocess, initargs=(stderr_lock, output_lock))
        with pool:
            tqdm_pathos.map(lambda task: design_task(*task), tasks, pool=pool)
    else:
        import tqdm
        from contextlib import nullcontext
        global STDERR_LOCK
        STDERR_LOCK = nullcontext()
        global OUTPUT_LOCK
        OUTPUT_LOCK = nullcontext()
        for task in tqdm.tqdm(tasks, desc="designing sequences"):
            design_task(*task)


def main():
    args = parse_args()
    from idrfeatlib import FeatureVector
    from idrfeatlib.utils import read_nested_csv, iter_nested, read_fasta, read_regions_csv
    from idrfeatlib.featurizer import compile_featurizer
    from idrfeatlib.native import compile_native_featurizer
    from idrfeatlib.metric import Metric
    from idrfeatlib.designer import FeatureDesigner, GreedyFeatureDesigner
    import os
    import json
    import sys
    import random
    import csv
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
    LENGTH_THRESHOLD = 30
    SEED_COLNAME = "Seed"
    MAX_SEED = 2 ** 64
    CONVERGENCE_THRESHOLD = 1e-4
    GOOD_MOVES_THRESHOLD = 3
    DECENT_MOVES_THRESHOLD = 5
    QUERY_COLNAME = "Sequence"
    if args.greedy:
        designer = GreedyFeatureDesigner(featurizer, metric, convergence_threshold=CONVERGENCE_THRESHOLD)
        designer.rng = random.Random() # type: ignore
    else:
        designer = FeatureDesigner(featurizer, metric, covergence_threshold=CONVERGENCE_THRESHOLD, good_moves_threshold=GOOD_MOVES_THRESHOLD, decent_moves_threshold=DECENT_MOVES_THRESHOLD, rng=random.Random())
    
    fa = dict(read_fasta(args.input_sequences))
    tasks = []
    colnames = ["ProteinID"]
    featnames = featurizer.keys()
    
    if args.input_regions is None:
        colnames += ["DesignID", "Sequence", "Time"]
        if args.save_seed:
            colnames.append("Seed")
        if args.keep_trajectory:
            colnames.append("Iteration")
        colnames += featnames
        fa = {protid: seq for protid, seq in fa.items() if len(seq) >= LENGTH_THRESHOLD}
        if args.query_file is None:
            if args.seeds_file is None:
                n_random = args.n_random or 1
                rng = random.Random()
                seeds = {protid: [rng.randint(0, MAX_SEED) for _ in range(n_random)] for protid in fa.keys()}
            else:
                seeds = read_nested_csv(args.seeds_file, 1, group_multiple=True)
                seeds = {protid: [row[SEED_COLNAME] for row in rows] for protid, rows in seeds.items() if protid in fa}
            for protid, prot_seeds in seeds.items():
                if (entry := fa.get(protid)) is None:
                    continue
                target = entry
                assert isinstance(target, str)
                for counter, seed in enumerate(prot_seeds):
                    seed = int(seed)
                    design_id = args.design_id.format(counter=counter, seed=seed, proteinid=protid)
                    tasks.append(
                        (None, target, protid, None, design_id, seed, designer, colnames, args)
                    )
        else:
            qries = read_nested_csv(args.query_file, 2)
            for protid, designid, row in iter_nested(qries, 2):
                if (entry := fa.get(protid)) is None:
                    continue
                target = entry
                assert isinstance(target, str)
                query = row[QUERY_COLNAME]
                tasks.append(
                    (query, target, protid, None, designid, None, designer, colnames, args)
                )
    else:
        colnames += ["RegionID", "DesignID", "Sequence", "Time"]
        if args.save_seed:
            colnames.append("Seed")
        if args.keep_trajectory:
            colnames.append("Iteration")
        colnames += featnames
        regions = read_regions_csv(args.input_regions)
        regions = {protid: ret for protid, entry in regions.items() if (ret := {regionid: (start, stop) for regionid, (start, stop) in entry.items() if stop - start >= LENGTH_THRESHOLD})}
        if args.query_file is None:
            if args.seeds_file is None:
                n_random = args.n_random or 1
                rng = random.Random()
                seeds = {protid: {regionid: [rng.randint(0, MAX_SEED) for _ in range(n_random)] for regionid in entry.keys()} for protid, entry in regions.items() if protid in fa}
                
            else:
                seeds = read_nested_csv(args.seeds_file, 2, group_multiple=True)
                seeds = {protid: {regionid: [row[SEED_COLNAME] for row in rows] for regionid, rows in entry.items()} for protid, entry in seeds.items() if protid in fa}
            for protid, regionid, region_seeds in iter_nested(seeds, 2):
                if (entry := fa.get(protid)) is None:
                    continue
                target_whole = entry
                assert isinstance(target_whole, str)
                start, stop = regions[protid][regionid] 
                target = target_whole[start:stop]
                if len(target) != stop - start:
                    print("invalid region `%s` for protein `%s` (start=%d,stop=%d,seqlen=%s)" % (regionid, protid, start, stop, len(target_whole)))
                    continue
                for counter, seed in enumerate(region_seeds):
                    seed = int(seed)
                    design_id = args.design_id.format(counter=counter, seed=seed, proteinid=protid, regionid=regionid, start=start, stop=stop)
                    tasks.append(
                        (None, target, protid, regionid, design_id, seed, designer, colnames, args)
                    )
        else:
            qries = read_nested_csv(args.query_file, 3)
            for protid, regionid, designid, row in iter_nested(qries, 3):
                if (entry := fa.get(protid)) is None:
                    continue
                target_whole = entry
                assert isinstance(target_whole, str)
                start, stop = regions[protid][regionid] 
                target = target_whole[start:stop]
                if len(target) != stop - start:
                    print("invalid region `%s` for protein `%s` (start=%d,stop=%d,seqlen=%s)" % (regionid, protid, start, stop, len(target_whole)))
                    continue
                query = row[QUERY_COLNAME]
                tasks.append(
                    (query, target, protid, regionid, designid, None, designer, colnames, args)
                )
    if not os.path.exists(args.output_file):
        with open(args.output_file, "w") as file:
            csv.DictWriter(file, colnames).writeheader()
        design_all(args.n_processes, tasks)
        return
    with open(args.output_file, "r") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames is None:
            with open(args.output_file, "w") as file:
                csv.DictWriter(file, colnames).writeheader()
            design_all(args.n_processes, tasks)
            return
        if reader.fieldnames != colnames:
            BOLD_RED = "\033[1;31m"
            NORMAL = "\033[0m"
            print(BOLD_RED + "cannot overwrite file `%s` with different column names (shown below):" % args.output_file + NORMAL, file=sys.stderr)
            print(",".join(colnames), file=sys.stderr)
            sys.exit(1)
        if args.keep_trajectory:
            checkpoint = [
                row for row in reader if row.pop("Iteration") == "END"
            ]
        else:
            checkpoint = list(reader)
    checkpoint_keys = ["ProteinID", "DesignID"]
    keys = [2, 4]
    if args.input_regions is not None:
        checkpoint_keys.insert(1, "RegionID")
        keys.insert(1, 3)
    checkpoint = {
        tuple(row.pop(key) for key in checkpoint_keys): row for row in checkpoint
    }
    tasks_not_done = []
    for task in tasks:
        checkpoint_key = tuple(task[key] for key in keys)
        if checkpoint_key not in checkpoint:
            tasks_not_done.append(task)
    design_all(args.n_processes, tasks_not_done)    


if __name__ == "__main__":
    main()