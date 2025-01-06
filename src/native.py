"""Module of default/native feature functions."""
import re
import typing
__all__ = [
    "compile_native_feature",
    "score_pattern_matches",
    "count_pattern_matches",
    "pattern_match_span",
    "log_ratio",
    "scd",
    "simple_spacing_closure",
    "custom_kappa_closure",
    "complexity",
    "isoelectric_point"
]
def compile_native_feature(
    *, 
    compute: str, 
    residue_groups: typing.Dict[str, str], 
    **kwargs
    ):
    """
    Turn kwargs into a sequence to feature function.
    
    You can use the kwargs directly. More often than not however, you'll use `**config`.
    """
    from functools import partial

    if compute == "score":
        if (score := kwargs.get("score")) is None:
            raise ValueError("`compute=score` requires a `score` parameter (see `score_pattern_matches`)")
        average = kwargs.get("take_average") or kwargs.get("average") or False
        return partial(score_pattern_matches, score=score, average=average)
    
    if compute == "count":
        if (pattern := kwargs.get("pattern")) is None:
            raise ValueError("`compute=count` requires a `pattern` parameter (see `count_pattern_matches`)")
        average = kwargs.get("take_average") or kwargs.get("average") or False
        return partial(count_pattern_matches, pattern=re.compile(pattern), average=average)
    
    if compute == "percent_residue":
        if (residue := kwargs.get("residue")) is None:
            raise ValueError("`compute=percent_residue` requires a residue as the `residue` parameter")
        return partial(count_pattern_matches, pattern=re.compile(residue), average=True)
    
    if compute == "percent_res_group":
        if (res_group_name := kwargs.get("residue_group")) is None:
            raise ValueError("`compute=percent_res_group` requires the name of a residue group as the `residue_group` parameter")
        if (res_group := residue_groups.get(res_group_name)) is None:
            raise ValueError("unknown residue group %s - available are: %s" % (res_group_name, ",".join(residue_groups.keys())))
        return partial(count_pattern_matches, pattern=re.compile("[%s]" % res_group), average=True)
    
    if compute == "span":
        if (pattern := kwargs.get("pattern")) is None:
            raise ValueError("`compute=span` requires a `pattern` parameter (see `pattern_match_span`)")
        return partial(pattern_match_span, pattern=re.compile(pattern))

    if compute == "log_ratio":
        if (num_aa := kwargs.get("numerator")) is None:
            raise ValueError("`compute=log_ratio` requires a `numerator` parameter (see `log_ratio`)") 
        if (denom_aa := kwargs.get("denominator")) is None:
            raise ValueError("`compute=log_ratio` requires a `denominator` parameter (see `log_ratio`)")
        return partial(log_ratio, num_aa=num_aa, denom_aa=denom_aa)
    
    if compute == "scd":
        return scd
    
    if compute == "simple_spacing":
        if (res_group_name := kwargs.get("residue_group")) is None:
            raise ValueError("`compute=simple_spacing` requires the name of a residue group as the `residue_group` parameter")
        if (res_group := residue_groups.get(res_group_name)) is None:
            raise ValueError("unknown residue group %s - available are: %s" % (res_group_name, ",".join(residue_groups.keys())))
        return simple_spacing_closure(res_group, res_group_name)
    
    if compute == "custom_kappa":
        return custom_kappa_closure()

    if compute == "complexity":
        return complexity
    
    if compute == "isoelectric_point":
        return isoelectric_point
    
    raise ValueError("not a recognized compute option: %s" % compute)

def score_pattern_matches(
    sequence: str,
    score: dict[re.Pattern[str], float],
    average: bool = False,
) -> float:
    """Calculate a weighted count or average of regex occurrences.

    Parameters
    ----------
    sequence : str
        Target sequence on which to perform the weighted count of
        regex occurrences.

    scores : dict[str, float]
        A dictionary containing the regex patterns to look for and the weights
        they contribute to the count.

    average : bool, optional
        Whether to divide by sequence length at the end.

        Defaults to ``False``.

    Raises
    ------
    If `average` is ``True`` and the provided sequence is empty.
    """
    result = sum(
        score * len(re.findall(pat, sequence)) for pat, score in score.items()
    )
    if average:
        return result / len(sequence)
    return result

def count_pattern_matches(
    sequence: str,
    pattern: re.Pattern[str],
    average: bool = False,
) -> float:
    """Count or average the number of regex occurrences in a sequence.

    Parameters
    ----------
    sequence : str
        Target sequence on which to count the regex occurrences.

    pattern : str
        The regex pattern to count.

    average : bool, optional
        Whether to divide by sequence length at the end.

        Defaults to ``False``.

    Raises
    ------
    If `average` is ``True`` and the provided sequence is empty.
    """
    result = len(re.findall(pattern, sequence))
    if average:
        return result / len(sequence)
    return result

def pattern_match_span(
    sequence: str,
    pattern: re.Pattern[str],
) -> float:
    """Calculate the total length spanned by patterns in a target sequence.

    Parameters
    ----------
    sequence : str
        Target sequence on which to determine the spanning length of the regex
        occurrences.

    pattern : str
        The regex pattern to determine the length of.
    """
    return sum(
        right - left for right, left in map(re.Match.span, re.finditer(pattern, sequence))
    )

def log_ratio(
    sequence: str,
    num_aa: str,
    denom_aa: str,
) -> float:
    from math import log
    """Calculate ``log(1 + num_aa) - log(1 - denom_aa)`` for in a sequence.

    Uses natural log.

    Parameters
    ----------
    sequence : str
        Target sequence on which to determine the log ratio.

    num_aa : str
        The amino acid in the numerator.

    denom_aa : str
        The amino acid in the denominator.
    """
    return log((1 + sequence.count(num_aa)) / (1 + sequence.count(denom_aa)))

def scd(sequence: str) -> float:
    """
    Calculate the `SCD`_ (Sequence Charge Decoration) of a sequence.

    .. _SCD: https://doi.org/10.1063/1.5005821

    Parameters
    ----------
    sequence : str
        Target sequence on which to determine the SCD.
    """
    from math import sqrt
    BINARY_CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
    charged_res = []
    for i, aa in enumerate(sequence):
        if aa in BINARY_CHARGE:
            charged_res.append(i)
    if not charged_res:
        return 0
    result: float = 0
    for i, loc_i in enumerate(charged_res):
        for loc_j in charged_res[:i]:
            result += (
                BINARY_CHARGE[sequence[loc_i]]
                * BINARY_CHARGE[sequence[loc_j]]
                * sqrt(loc_i - loc_j)
            )
    return result / len(sequence)

def abstract_spacing_calculation(
    *,
    sequence: str,
    candidates: typing.List[int],
    are_neighbours: typing.Callable[[str, str], bool],
    prob_neighbor_given_candidate: typing.Callable[
        [str, typing.List[int], int], float
    ],
    not_enough_candidates_error: str,
    blob: int,
) -> float:
    """Calculate a spacing parameter.

    This is an abstraction over two similar computations: `simple_spacing`
    and `custom_kappa`. They both involve counting two subsets of residues:
    1. `candidates`, given by the provided list, and
    2. `neighbours`, which are a subset of those candidates that must be
        within close proximity to each other.

    The resulting number is a measure of how often close-neighbour residue
    pairs occur compared to how often they are expected to occur via
    composition.
            
    This function is not publically exposed as a `*` import.
    """
    from math import sqrt
    def count_neighbors():
        """The number of `neighbour` residues."""
        candidate_pairs = list(zip(candidates[:-1], candidates[1:]))
        if not candidate_pairs:
            raise ValueError(not_enough_candidates_error)
        return sum(
            1
            for i, j in candidate_pairs
            if are_neighbours(sequence[i], sequence[j]) and (abs(i - j) <= blob)
        )
    actual_neighbors = count_neighbors()
    p: float = prob_neighbor_given_candidate(
        sequence, candidates, blob
    )
    mean_neighbors: float = p * len(candidates)
    sd_neighbors: float = sqrt(p * (1 - p) * len(candidates))
    return (actual_neighbors - mean_neighbors) / sd_neighbors


def simple_spacing_closure(
    res_group: str,
    res_group_name: str,
    *,
    blob: int = 5,
):
    """
    Set up a closure to compute simple spacing for this residue group.
    
    In a simple spacing computation, the residues in a certain residue group (e.g.
    the charged residues, the aromatic residues, etc.) are considered to be
    `candidate` residues. The pairs of candidate residues within 5 residues of each
    other are the `neighbour` residues.

    This closure therefore measures how "clustered" the residues from this residue
    group are distributed in the sequence.
    """
    def prob_neighbor_given_candidate(
        sequence: str,
        candidates: typing.List[int],
        blob: int,
    ):
        proportion_candidates: float = len(candidates) / len(sequence)
        # It shouldn't be zero either, but that should be caught already.
        if proportion_candidates == 1:
            raise ValueError("cannot compute %s spacing on sequence with only %s residues" % (res_group_name, res_group_name))
        return proportion_candidates * sum(
            (1 - proportion_candidates) ** i for i in range(blob)
        )
    def are_neighbours(_a: str, _b: str):
        return True
    not_enough_candidates_error = "sequence has no %s residues" % res_group_name
    def simple_spacing(sequence: str):
        candidates = []
        for i, aa in enumerate(sequence):
            if aa in res_group:
                candidates.append(i)
        return abstract_spacing_calculation(
            sequence=sequence,
            candidates=candidates,
            are_neighbours=are_neighbours,
            prob_neighbor_given_candidate=prob_neighbor_given_candidate,
            not_enough_candidates_error=not_enough_candidates_error,
            blob=blob
        )
    return simple_spacing


def custom_kappa_closure(
    *,
    blob: int = 5,
):
    """
    Set up a closure to compute a kappa-like measure for this residue group.
    
    In this custom-kappa computation, charged residues are candidates, and the pairs of
    same-charge residues within 5 residues of each other are the `neighbour` residues.

    This closure therefore measures how similar/blocky charges residues are.
    """
    BINARY_CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1} 
    def prob_neighbor_given_candidate(
        sequence: str,
        candidates: typing.List[int],
        blob: int,
    ):
        proportion_candidates: float = len(candidates) / len(sequence)
        count_pos: int = sum(1 for i in candidates if BINARY_CHARGE[sequence[i]] == 1)
        count_neg: int = len(candidates) - count_pos
        prob_charges_are_diff: float = (
            2 * (count_pos) * (count_neg) / (len(candidates) ** 2)
        )
        if proportion_candidates == 1 and prob_charges_are_diff == 0:
            raise ValueError("cannot compute charge spacing, `custom_kappa` cannot deal with all residues of one charge")
        prob_next_charge_in_blob: float = proportion_candidates * sum(
            (1 - proportion_candidates) ** i for i in range(blob)
        )
        return prob_next_charge_in_blob * (1 - prob_charges_are_diff)
    def are_neighbours(res_a: str, res_b: str):
        return BINARY_CHARGE[res_a] == BINARY_CHARGE[res_b]
    not_enough_candidates_error = "sequence has no charged residues"
    def custom_kappa(sequence: str):
        candidates = []
        for i, aa in enumerate(sequence):
            if aa in BINARY_CHARGE:
                candidates.append(i)
        return abstract_spacing_calculation(
            sequence=sequence,
            candidates=candidates,
            are_neighbours=are_neighbours,
            prob_neighbor_given_candidate=prob_neighbor_given_candidate,
            not_enough_candidates_error=not_enough_candidates_error,
            blob=blob
        ) 
    return custom_kappa

def complexity(sequence: str) -> float:
    """Calculate the `complexity`_ (entropy-like feature) of a target sequence.

    Uses natural log.

    .. _complexity: https://www.sciencedirect.com/science/article/pii/009784859385006X

    Parameters
    ----------
    sequence : str
        Target sequence on which to determine the complexity.

    Raises
    ------
    If the sequence is empty.
    """  # noqa: E501, pylint: disable=line-too-long
    from math import lgamma
    AMINOACIDS = "ACDEFGHIKLMNPQRSTVWY"
    log_gamma_sum: float = 0
    for aa in AMINOACIDS:
        log_gamma_sum += lgamma(1 + sequence.count(aa))
    return (lgamma(1 + len(sequence)) - log_gamma_sum) / len(sequence)

def accurate_net_charge(
    ph: float,
    num_basic_res: int,
    counts_and_pkas: typing.Iterable[tuple[int, float]],
) -> float:
    """Calculate a very accurate net charge based on pKa formulas.

    Helper for :func:`~isoelectric_point`.

    Parameters
    ----------
    ph : float
        The pH to be calculated at.

    num_basic_res: int
        The number of sites which are positively charged when protonated,
        including the basic N terminus site. i.e.
        0 is not a valid value because there is always a basic N terminus site.

    counts_and_pkas : Iterable[tuple[int, float]]
        The number of sites of a specific pKa, for each pKa.

    Raises
    ------
    If an invalid (non-positive) `num_basic_res` is passed in.

    Note
    ----
    **Math logic**

    The logic behind this function is that one counts the basic sites
    (positively charged in their protonated state) as the default charge.

    As you increase the pH, a proportion of those sites become deprotonated,
    decreasing the charge. The proportion of deprotonated sites of species
    with a fixed pka is:

    >>> 1 / (1 + 10 ** (pka - ph))

    So you subtract the expected number of free protons (corresponding to
    deprotonated sites) from the number of basic sites.
    """

    PKA_N_TERM = 7.5
    PKA_C_TERM = 3.55

    if num_basic_res < 1:
        raise AssertionError(
            "Trying to calculate net charge on "
            + "a negative number of basic residues! (%d)" % num_basic_res,
        )
    free_protons: float = 0
    for count, pka in counts_and_pkas:
        proportion_protonated = 1 / (1 + 10 ** (ph - pka))
        free_protons += count * (1 - proportion_protonated)
    proportion_protonated = 1 / (1 + (10 ** (ph - PKA_N_TERM)))
    free_protons += 1 - proportion_protonated
    proportion_protonated = 1 / (1 + (10 ** (ph - PKA_C_TERM)))
    free_protons += 1 - proportion_protonated
    return num_basic_res - free_protons


def binary_search_root_finder(
    f: typing.Callable[[float], float],
    bracket: tuple[float, float],
    threshold: float = 10 * -4,
) -> float:
    """Find roots of a decreasing sigmoidal function.

    Dependency for `isoelectric_point`.

    Parameters
    ----------
    f : Callable[[float], float]
        Function to find the root of. Should only have one root.

    guess : float, optional
        Initial guess.

    bracket : tuple[float, float]
        Interval on which to find the root.
    """
    bottom, top = bracket
    guess = (top + bottom) / 2
    while top - bottom > threshold:
        if f(guess) > 0:
            bottom = guess
        else:
            top = guess
        guess = (top + bottom) / 2
    return guess

def isoelectric_point(sequence: str):
    """Calculate the isoelectric point of a target sequence.

    Searches for a root of the pH-charge curve on the pH interval 0 to 14.

    Parameters
    ----------
    sequence : str
        Target sequence on which to determine the isoelectric point.

    Raises
    ------
    If the charge curve is negative at the acidic end or positive at the
    basic end, it is guaranteed not to have a root on the interval.
    """
    from functools import partial

    PKAS_ALL = {
        "D": 4.05,
        "E": 4.45,
        "C": 9.0,
        "Y": 10.0,
        "K": 10.0,
        "R": 12.0,
        "H": 5.98
    }
    BASIC_RES = "KHR"

    counts_and_pkas: list[tuple[int, float]] = [
        (sequence.count(aa), pka) for aa, pka in PKAS_ALL.items()
    ]
    num_basic_res: int = 1 + sum(1 for aa in sequence if aa in BASIC_RES)

    continuous_charge: typing.Callable[[float], float] = partial(
        accurate_net_charge,
        num_basic_res=num_basic_res,
        counts_and_pkas=counts_and_pkas,
    )
    if (continuous_charge(0) < 0) or (continuous_charge(14) > 0):
        raise ValueError("Isoelectric point of the protein "
            + "cannot be on the interval [0, 14].")
    return binary_search_root_finder(continuous_charge, (0, 14))