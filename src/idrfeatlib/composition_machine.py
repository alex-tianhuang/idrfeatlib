"""
Module containing the CompositionMachine class,
a class that optimizes amino acid composition over common functions.

Uses numpy.
"""
import typing
__all__ = [
    "CompositionMachine"
]

class CompositionMachine:
    def __init__(self,
        sum_scores: typing.Dict[str, typing.Dict[str, float]],
        sum_score_weights: typing.Dict[str, float],
        sum_score_targets: typing.Dict[str, float],
        average_scores: typing.Dict[str, typing.Dict[str, float]],
        average_score_weights: typing.Dict[str, float],
        average_score_targets: typing.Dict[str, float],
        percent_resgroup_weights: typing.Dict[str, float],
        percent_resgroup_targets: typing.Dict[str, float],
        log_ratio_weights: typing.Dict[typing.Tuple[str, str], float],
        log_ratio_targets: typing.Dict[typing.Tuple[str, str], float],
        complexity_weight: float,
        complexity_target: float,
        isoelectric_point_weight: float,
        isoelectric_point_target: float,
    ) -> None:
        self.sum_scores = sum_scores
        self.sum_score_weights = sum_score_weights
        self.sum_score_targets = sum_score_targets
        self.average_scores = average_scores
        self.average_score_weights = average_score_weights
        self.average_score_targets = average_score_targets
        self.percent_resgroup_weights = percent_resgroup_weights
        self.percent_resgroup_targets = percent_resgroup_targets
        self.log_ratio_weights = log_ratio_weights
        self.log_ratio_targets = log_ratio_targets
        self.complexity_weight = complexity_weight
        self.complexity_target = complexity_target
        self.isoelectric_point_weight = isoelectric_point_weight
        self.isoelectric_point_target = isoelectric_point_target

    def autograd_machine(self, length: int):
        """
        Make a closure that simultaneously computes the weighted squared deviation and amino acid gradients.
        
        Sequence complexity gradient is approximated.

        Dev note
        --------
        The structure of this function is a little unorthodox because typically Python people like `_private` functions
        But I prefer local closures defined in the body of the function because underscores are ugly.
        
        Most of the function is defining various gradient computing closures for each subtype of feature
        and then bundling them all together in an `autograd_closure_list`, which acts like a virtualized
        feature function vector.
        """
        import math
        def autograd_score_closure(scoring: typing.Dict[str, float]):
            gradient = scoring
            def autograd_score(composition: typing.Dict[str, int]):
                nonlocal gradient
                score = sum(
                    scoring.get(aa, 0) * count for aa, count in composition.items()
                )
                return score, gradient
            return autograd_score
        score_closures = [(autograd_score_closure(scoring), self.sum_score_targets[score_name], self.sum_score_weights[score_name]) for score_name, scoring in self.sum_scores.items()] \
        + [(autograd_score_closure({aa: score / length for aa, score in scoring.items()}), self.average_score_targets[score_name], self.average_score_weights[score_name]) for score_name, scoring in self.average_scores.items()] 
        def autograd_percent_resgroup_closure(res_group: str):
            gradient = {
                aa: 1 / length for aa in res_group
            }
            def autograd_percent_resgroup(composition: typing.Dict[str, int]):
                nonlocal gradient, length
                return sum(composition[aa] for aa in res_group) / length, gradient
            return autograd_percent_resgroup
        percent_resgroup_closures = [(autograd_percent_resgroup_closure(res_group), target, self.percent_resgroup_weights[res_group]) for res_group, target in self.percent_resgroup_targets.items()]
        def autograd_log_ratio_closure(num_aa: str, denom_aa: str):
            def autograd_log_ratio(composition: typing.Dict[str, int]):
                numerator = 1 + composition[num_aa]
                denominator = 1 + composition[denom_aa]
                gradient = {
                    num_aa: 1 / numerator,
                    denom_aa: - 1 / denominator
                }
                return math.log(numerator / denominator), gradient
            return autograd_log_ratio
        log_ratio_closures = [(autograd_log_ratio_closure(num_aa, denom_aa), target, self.log_ratio_weights[(num_aa, denom_aa)]) for (num_aa, denom_aa), target in self.log_ratio_targets.items()]
        log_factorial_length = math.lgamma(1 + length)
        def autograd_complexity(composition: typing.Dict[str, int]):
            nonlocal log_factorial_length, length
            log_gamma_sum: float = 0
            gradient = {}
            for aa, count in composition.items():
                log_gamma_sum += math.lgamma(1 + count)
                # could use digamma as d[lgamma]/dx but that seems expensive
                gradient[aa] = - math.log(1 + count) / length 
            return (log_factorial_length - log_gamma_sum) / length, gradient
        PKAS_ALL = {
            "D": 4.05,
            "E": 4.45,
            "C": 9.0,
            "Y": 10.0,
            "K": 10.0,
            "R": 12.0,
            "H": 5.98,
        }
        BASIC_RES = "KHR"
        ACIDIC_RES = "DECY"
        def gradients_accurate_net_charge(ph: float, composition: typing.Dict[str, int]):
            nonlocal ACIDIC_RES, PKAS_ALL
            aa_gradient = {
                **{aa: -1.0 for aa in ACIDIC_RES},
                **{aa: 0.0 for aa in BASIC_RES}
            }
            ph_derivative = 0
            for aa, pka in PKAS_ALL.items():
                proportion_protonated = 1 / (1 + 10 ** (ph - pka))
                aa_gradient[aa] += proportion_protonated
                ph_derivative -= composition[aa] * (proportion_protonated - proportion_protonated * proportion_protonated)
            ph_derivative *= math.log(10)
            return aa_gradient, ph_derivative
        def gradient_isoelectric_point(ph: float, composition: typing.Dict[str, int]):
            # d[CH] = d[CH]/d[AA] x d[AA] + d[CH]/d[PH] x d[PH]
            # 0 = d[CH]/d[AA] x d[AA] + d[CH]/d[PH] x d[PH]
            # d[CH]/d[AA] x d[AA] = - d[CH]/d[PH] x d[PH]
            # d[CH]/d[PH] x d[PH] = - d[CH]/d[AA] x d[AA]
            # d[PH]/d[AA] = - d[CH]/d[AA] / d[CH]/d[PH]
            aa_gradient, ph_derivative = gradients_accurate_net_charge(ph, composition)
            return {
                aa: - g / ph_derivative for aa, g in aa_gradient.items()
            }
        def autograd_isoelectric_point(composition: typing.Dict[str, int]):
            from .native import binary_search_root_finder, accurate_net_charge
            from functools import partial
            nonlocal PKAS_ALL, BASIC_RES, gradient_isoelectric_point
            counts_and_pkas: list[tuple[int, float]] = [
                (composition[aa], pka) for aa, pka in PKAS_ALL.items()
            ]
            num_basic_res: int = 1 + sum(composition[aa] for aa in BASIC_RES)

            continuous_charge: typing.Callable[[float], float] = partial(
                accurate_net_charge,
                num_basic_res=num_basic_res,
                counts_and_pkas=counts_and_pkas,
            )
            if (continuous_charge(0) < 0) or (continuous_charge(14) > 0):
                raise ValueError(
                    "Isoelectric point of the protein " + "cannot be on the interval [0, 14]."
                )
            ph = binary_search_root_finder(continuous_charge, (0, 14))
            gradient = gradient_isoelectric_point(ph, composition)
            return ph, gradient
        autograd_closure_list = score_closures + percent_resgroup_closures + log_ratio_closures
        if self.complexity_weight != 0.0:
            autograd_closure_list.append((autograd_complexity, self.complexity_target, self.complexity_weight))
        if self.isoelectric_point_weight != 0.0:
            autograd_closure_list.append((autograd_isoelectric_point, self.isoelectric_point_target, self.isoelectric_point_weight))
        def autograd(composition: typing.Dict[str, int]):
            nonlocal autograd_closure_list
            ssr = 0
            sum_gradient = {
                aa: 0.0 for aa in composition.keys()
            }
            for func, target, weight in autograd_closure_list:
                score, gradient = func(composition)
                residual = score - target
                chain_rule_factor = weight * residual
                ssr += chain_rule_factor * residual
                for aa, grad_value in gradient.items():
                    sum_gradient[aa] += chain_rule_factor * grad_value
            return ssr, sum_gradient
        return autograd

    def search_optimal_composition(self, length: int):
        """
        Search the composition space greedily for an optimal composition solution.
        """
        from itertools import product
        AMINOACIDS = "ACDEFGHIKLMNPQRSTVWY"
        remaining = length
        guess = {
            aa: 0 for aa in AMINOACIDS
        }
        chosen_count = round(length / 20)
        for aa in AMINOACIDS:
            if chosen_count > remaining:
                guess[aa] = remaining
                break
            guess[aa] = chosen_count
            remaining -= chosen_count
        if remaining > 0:
            guess[aa] += remaining
        
        autograd = self.autograd_machine(length)
        ssr, gradient = autograd(guess)
        NUM_TRIES = 100 * length
        for _ in range(NUM_TRIES):
            gradient_sorted_aas = sorted(AMINOACIDS, key=lambda aa: gradient[aa])
            dec_aas = [aa for aa in reversed(gradient_sorted_aas) if guess[aa] > 0]
            inc_aas = [aa for aa in gradient_sorted_aas if guess[aa] < length]
            best_move = None
            for dec_aa, inc_aa in product(dec_aas, inc_aas):
                if dec_aa == inc_aa:
                    continue
                guess[dec_aa] -= 1
                guess[inc_aa] += 1
                new_ssr, gradient = autograd(guess)
                if new_ssr < ssr:
                    ssr = new_ssr
                    best_move = (dec_aa, inc_aa)
                    break
                guess[dec_aa] += 1
                guess[inc_aa] -= 1 
            else:
                if best_move is None:
                    return guess
                dec_aa, inc_aa = best_move # type: ignore
                guess[dec_aa] -= 1
                guess[inc_aa] += 1
                ssr, gradient = autograd(guess)
            if ssr == 0.0:
                return guess
        raise RuntimeError("Composition search failed to converge after {} iterations.".format(NUM_TRIES))

    @staticmethod
    def new(weights: typing.Dict[str, float], targets: typing.Dict[str, float], features_dict: typing.Any = None):
        """
        Make a CompositionMachine from weights and targets.
        
        This is the preferred way to make a composition machine.
        """
        import json
        import os
        if features_dict is None:
            with open(
                os.path.join(os.path.dirname(__file__), "native-features.json"), "r"
            ) as file:
                features_dict = json.load(file)
        return_value = CompositionMachine(
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            0,
            0,
            0,
            0
        )
        def compile_one_feature(
                featname: str, 
                compute: str, 
                residue_groups: typing.Dict[str, str],
                # motif_frequencies: typing.Dict[str, float],
                # residue_frequencies: typing.Optional[typing.Dict[str, float]],
                **kwargs
            ):
            if compute == "score":
                if (score := kwargs.get("score")) is None:
                    raise ValueError(
                        "`compute=score` requires a `score` parameter (see `score_pattern_matches`)"
                    )
                average = kwargs.get("take_average") or kwargs.get("average") or False
                if not isinstance(average, bool):
                    raise TypeError("expected `average` to be True or False")
                if not isinstance(score, dict):
                    raise TypeError("expected `score` to be a dict of residue->score pairs")
                for key in score.keys():
                    if len(key) != 1:
                        raise ValueError("expected `score` to be a dict of residue->score pairs, got a non-residue string")
                if average:
                    return_value.average_scores[featname] = score
                    return_value.average_score_weights[featname] = weights[featname]
                    return_value.average_score_targets[featname] = targets[featname]
                else:
                    return_value.sum_scores[featname] = score
                    return_value.sum_score_weights[featname] = weights[featname]
                    return_value.sum_score_targets[featname] = targets[featname]
                return
            

            if compute == "percent_residue":
                if (residue := kwargs.get("residue")) is None:
                    raise ValueError(
                        "`compute=percent_residue` requires a residue as the `residue` parameter"
                    )
                if not isinstance(residue, str) and len(residue) == 1:
                    raise TypeError("expected `residue` to be a single amino acid")
                return_value.percent_resgroup_weights[residue] = weights[featname]
                return_value.percent_resgroup_targets[residue] = targets[featname]
                return

            if compute == "percent_res_group":
                if (res_group_name := kwargs.get("residue_group")) is None:
                    raise ValueError(
                        "`compute=percent_res_group` requires the name of a residue group as the `residue_group` parameter"
                    )
                if (res_group := residue_groups.get(res_group_name)) is None:
                    raise ValueError(
                        "unknown residue group %s - available are: %s"
                        % (res_group_name, ",".join(residue_groups.keys()))
                    )
                if isinstance(res_group, list):
                    res_group = "".join(res_group)
                if not isinstance(res_group, str):
                    raise TypeError(
                        "expected residue group %s to be a list or string of amino acids"
                        % res_group_name
                    )
                return_value.percent_resgroup_weights[res_group] = weights[featname]
                return_value.percent_resgroup_targets[res_group] = targets[featname]
                return

            if compute == "log_ratio":
                if (num_aa := kwargs.get("numerator")) is None:
                    raise ValueError(
                        "`compute=log_ratio` requires a `numerator` parameter (see `log_ratio`)"
                    )
                if (denom_aa := kwargs.get("denominator")) is None:
                    raise ValueError(
                        "`compute=log_ratio` requires a `denominator` parameter (see `log_ratio`)"
                    )
                return_value.log_ratio_weights[(num_aa, denom_aa)] = weights[featname]
                return_value.log_ratio_targets[(num_aa, denom_aa)] = targets[featname]
                return

            if compute == "complexity" or compute == "sequence_complexity":
                return_value.complexity_weight = weights[featname]
                return_value.complexity_target = targets[featname]
                return

            if compute == "isoelectric_point":
                return_value.isoelectric_point_weight = weights[featname]
                return_value.isoelectric_point_target = targets[featname]
                return

        features = features_dict["features"]
        residue_groups = features_dict.get("residue_groups") or {}
        # motif_frequencies = features_dict.get("motif_frequencies") or {}
        # residue_frequencies = features_dict.get("aa_frequencies")
        errors = {}
        for featname, feature_params in features.items():
            try:
                compile_one_feature(
                    featname=featname,
                    residue_groups=residue_groups,
                    # motif_frequencies=motif_frequencies,
                    # residue_frequencies=residue_frequencies,
                    **feature_params,
                )
            except (ValueError, TypeError) as e:
                errors[featname] = e
        return return_value, errors
       
        
        