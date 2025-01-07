"""
Contains the `Featurizer` class, a convenience class for running sequence feature computations.

To make one of these, you need a dictionary of function pointers. To get the default/native
sequence feature functions, see `compile_native_featurizer` in the `native` module in the
same folder.

Otherwise, feel free to add your own function pointers.
"""
from . import FeatureVector
import typing

__all__ = ["Featurizer"]

class Featurizer:
    """
    An object that holds a dict of (sequence -> float) functions, i.e. `features`.
    
    Call `Featurizer.featurize` to compute the feature vector of a single sequence.
    Call `Featurizer.featurize_to_matrices` to compute a lot of feature vectors.
    
    To make one of these, see `compile_native_featurizer` in the `native` module in the
    same folder.
    """
    def __init__(self, funcs: typing.Dict[str, typing.Callable[..., float]]) -> None:
        """`funcs` should be a dict of (sequence -> float) functions."""
        self._funcs = funcs
    def featurize(self, sequence: str, *, acceptable_errors = (ArithmeticError, ValueError, KeyError)) -> typing.Tuple["FeatureVector", typing.Dict[str, Exception]]:
        """Compute the feature vector of a single sequence, and also return its failed computations."""
        return_value = {}
        errors = {}
        for featname, func in self._funcs.items():
            try:
                return_value[featname] = func(sequence)
            except acceptable_errors as e:
                errors[featname] = e
        return FeatureVector(return_value), errors
    def featurize_to_matrices(self, sequences: typing.Iterable[typing.Tuple[str, str]], *, acceptable_errors = (ArithmeticError, ValueError, KeyError)) -> typing.Tuple[typing.Dict[str, "FeatureVector"], typing.Dict[str, typing.Dict[str, Exception]]]:
        """Compute the feature vector of many sequences, and also return their failed computations."""
        fvecs_all = {}
        errors_all = {}
        for label, sequence in sequences:
            fvec, errors = self.featurize(sequence, acceptable_errors=acceptable_errors)
            fvecs_all[label] = fvec
            if errors:
                errors_all[label] = errors
        return fvecs_all, errors_all
