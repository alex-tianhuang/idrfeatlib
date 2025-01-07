import typing
from .metric import Metric
from .native import scd
from time import time
from .featurizer import Featurizer
from random import Random
from copy import copy

__all__ = [
    "FeatureDesigner",   
]
class FeatureDesigner:
    """Feature mimic design algorithm."""
    def __init__(self, featurizer: typing.Dict[str, typing.Callable[..., float]], metric: Metric, covergence_threshold: float, good_moves_threshold: int, rng: Random) -> None:
        self.featurizer = featurizer
        self.metric = metric
        self.convergence_threshold = covergence_threshold
        self.good_moves_threshold = good_moves_threshold
        self.rng = rng

    def generate_mutations(self, length: int):
        """Helper for `design_loop`."""
        AMINOACIDS = list("ACDEFGHIKLMNPQRSTVWY")
        N_AA = len(AMINOACIDS)
        seeds = list(range(length * N_AA))
        self.rng.shuffle(seeds)
        for seed in seeds:
            aa = AMINOACIDS[seed % N_AA]
            pos = int(seed / N_AA)
            yield pos, aa
        
    def design_loop(
        self,
        query: str,
        acceptable_errors = (ArithmeticError, ValueError, KeyError) 
    ):
        """
        Optimize a query sequence to fit the feature vector at `self.metric.origin`.

        Parameters
        ----------
        query : str
            The starting sequence to generate designs off of by iterative mutation.
        """
        query_fvec, _ = Featurizer(self.featurizer).featurize(query, acceptable_errors=())

        scd_machine = None
        scd_key = None
        featurizer = self.featurizer
        
        for featname, feature in list(self.featurizer.items()):
            if feature is scd:
                if scd_machine is None:
                    featurizer = copy(self.featurizer)
                    previous_scd = query_fvec.as_dict[featname]
                    charged_res = set()
                    for i, aa in enumerate(query):
                        if aa in CHARGE:
                            charged_res.add(i)
                    mutation_from = "" # placeholder for string type
                    mutation_pos = len(query) # should throw error if used
                    scd_machine = ScdMachine(previous_scd, charged_res, mutation_from, mutation_pos)
                    scd_key = featname
                featurizer[featname] = scd_machine.compute_scd

        featurizer = Featurizer(featurizer)

        curr_dist_to_tgt: float = self.metric.euclidean_norm_of(query_fvec)
        start_time = time()
        iteration = 0
        yield {
            **query_fvec.as_dict,
            "Iteration": iteration,
            "Sequence": query,
            "Time": 0
        }
        next_dist_to_tgt: float = curr_dist_to_tgt
        curr_dist_to_tgt += self.convergence_threshold + 1

        while curr_dist_to_tgt - next_dist_to_tgt > self.convergence_threshold:
            curr_dist_to_tgt = next_dist_to_tgt
            all_candidates = [(query, copy(scd_machine))]

            for mutation_pos, mutation_to in self.find_lower_mutations(
                self.generate_mutations(len(query)), curr_dist_to_tgt, query, featurizer, scd_machine, acceptable_errors=acceptable_errors
            ):
                new_candidates = []
                for old_candidate, old_scd_machine in all_candidates:
                    if old_scd_machine is not None:
                        assert scd_machine is not None
                        scd_machine.clone_from(old_scd_machine)
                    assert (candidate := apply_mutation(old_candidate,mutation_pos, mutation_to, scd_machine)) is not None
                    try:
                        candidate_fvec, _ = featurizer.featurize(candidate, acceptable_errors=())    
                    except acceptable_errors:
                        continue
                    candidate_dist_to_tgt = self.metric.euclidean_norm_of(candidate_fvec)
                    if next_dist_to_tgt > candidate_dist_to_tgt:
                        next_dist_to_tgt = candidate_dist_to_tgt
                        query_fvec = candidate_fvec
                        query = candidate
                    if scd_machine is not None:
                        assert scd_key is not None
                        new_scd = candidate_fvec.as_dict[scd_key]
                        scd_machine.advance_mutation(mutation_pos, mutation_to, new_scd)
                        new_scd_machine = scd_machine.shallow_clone()
                    else:
                        new_scd_machine = None
                    new_candidates.append((candidate, new_scd_machine))
                all_candidates.extend(new_candidates)
                
            yield {
                **query_fvec.as_dict,
                "Iteration": iteration,
                "Sequence": query,
                "Time": time() - start_time
            }

        yield {
            **query_fvec.as_dict,
            "Iteration": "END",
            "Sequence": query,
            "Time": 0
        }
        return query

    def find_lower_mutations(
        self,
        mutation_generator: typing.Iterator[typing.Tuple[int, str]],
        curr_dist_to_tgt: float,
        current_seq: str,
        featurizer: Featurizer,
        scd_machine: typing.Optional["ScdMachine"],
        acceptable_errors: ...,
    ) -> typing.List[typing.Tuple[int, str]]:
        """
        Find single point mutations that decrease the distance to a target.

        Bins them into those that decrease the distance to target a lot
        ("good_moves") and those that decrease the distance by a small amount
        ("decent_moves").
        """
        def add_mutation_to_memo(
            mutation_pos: int,
            mutation_to: str,
            dist_to_tgt: float,
            memo: dict[int, typing.Tuple[str, float]],
        ):
            """Helper for `find_lower_mutations`."""
            if (collision := memo.get(mutation_pos)) is None:
                memo[mutation_pos] = mutation_to, dist_to_tgt
                return
            _, memoized_dist = collision
            if dist_to_tgt >= memoized_dist:
                return
            memo[mutation_pos] = mutation_to, dist_to_tgt
        good_mutation_dists = {}
        decent_mutation_dists = {}

        for mutation_pos, mutation_to in mutation_generator:
            if (guess := apply_mutation(current_seq, mutation_pos, mutation_to, scd_machine)) is None:
                continue
            try:
                guess_fvec, _ = featurizer.featurize(guess, acceptable_errors=())
            except acceptable_errors:
                continue
            guess_dist_to_tgt = self.metric.euclidean_norm_of(guess_fvec)
            if guess_dist_to_tgt < curr_dist_to_tgt:
                if (
                    curr_dist_to_tgt - guess_dist_to_tgt
                    > self.convergence_threshold
                ):
                    add_mutation_to_memo(
                        mutation_pos, mutation_to, guess_dist_to_tgt, good_mutation_dists
                    )
                else:
                    add_mutation_to_memo(
                        mutation_pos, mutation_to, guess_dist_to_tgt, decent_mutation_dists
                    )
            if len(good_mutation_dists) >= self.good_moves_threshold:
                break
        mutations = {
            pos: aa for pos, (aa, _) in decent_mutation_dists.items()
        } 
        mutations.update({
            pos: aa for pos, (aa, _) in good_mutation_dists.items() 
        })
        return list(mutations.items())

CHARGE = {
    "D": -1,
    "E": -1,
    "K": 1,
    "R": 1
}

def apply_mutation(sequence: str, mutation_pos: int, mutation_to: str, scd_machine: typing.Optional["ScdMachine"]):
    """Apply a point mutation to a sequence, appropriately adjusting the `scd_machine`."""
    if (mutation_from := sequence[mutation_pos]) == mutation_to:
        return None
    if scd_machine is not None:
        scd_machine.mock_mutation(mutation_from, mutation_pos)
    return sequence[:mutation_pos] + mutation_to + sequence[mutation_pos+1:]
    
class ScdMachine:
    """A helper for the mimic design algorithm that computes SCD in non-quadratic time."""
    def __init__(self, previous_scd: float, charged_res: typing.Set[int], mutation_from: str, mutation_pos: int) -> None:
        self.previous_scd = previous_scd
        self.charged_res = charged_res
        self.mutation_from = mutation_from
        self.mutation_pos = mutation_pos
    def compute_scd(self, sequence: str):
        """Compute the SCD from a previous SCD and mutation data."""
        from math import sqrt
        ch_delta = CHARGE.get(sequence[self.mutation_pos], 0) - CHARGE.get(self.mutation_from, 0)
        if ch_delta == 0:
            return self.previous_scd
        scd_delta = 0
        for i in self.charged_res:
            charge_at_i = CHARGE[sequence[i]]
            scd_delta = ch_delta * charge_at_i * sqrt(abs(self.mutation_pos - i))
        return self.previous_scd + scd_delta / len(sequence)
    def mock_mutation(self, mutation_from: str, mutation_pos: int):
        """
        Update the mutation data, but not the SCD data.
        
        This is called before a `featurizer.featurize` call.
        """
        self.mutation_from = mutation_from
        self.mutation_pos = mutation_pos
    def advance_mutation(self, mutation_pos: int, mutation_to: str, scd: float):
        """
        Update the SCD data.

        This is called before the next round of `mock_mutation`s.
        """
        if mutation_to not in CHARGE:
            self.charged_res.remove(mutation_pos)
        else:
            self.charged_res.add(mutation_pos)
        self.previous_scd = scd
    def shallow_clone(self):
        return ScdMachine(self.previous_scd, self.charged_res, self.mutation_from, self.mutation_pos)
    def clone(self):
        """Deepcopy."""
        return ScdMachine(self.previous_scd, copy(self.charged_res), self.mutation_from, self.mutation_pos)
    def clone_from(self, other):
        """Deepcopy from other into self."""
        self.previous_scd = other.previous_scd
        self.charged_res = copy(other.charged_res)
        self.mutation_from = other.mutation_from
        self.mutation_pos = other.mutation_pos
