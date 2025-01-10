"""
Contains the `Featurizer` class, a convenience class for running sequence feature computations.

To make one of these, you need a dictionary of function pointers. To get the default/native
sequence feature functions, see `compile_native_featurizer` in the `native` module in the
same folder.

Also feel free to add your own function pointers either using the `compile_featurizer` in this
module, or just by constructing your own dict of named function pointers.
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
    def featurize_to_matrices(self, sequences: typing.Iterable[typing.Tuple[typing.Any, str]], *, acceptable_errors = (ArithmeticError, ValueError, KeyError)) -> typing.Tuple[typing.Dict[typing.Any, "FeatureVector"], typing.Dict[typing.Any, typing.Dict[str, Exception]]]:
        """
        Compute the feature vector of many sequences, and also return their failed computations.
        
        The two returned dicts (feature vector and errors) are indexed by whatever `sequences` was indexed by.
        """
        fvecs_all = {}
        errors_all = {}
        for label, sequence in sequences:
            fvec, errors = self.featurize(sequence, acceptable_errors=acceptable_errors)
            fvecs_all[label] = fvec
            if errors:
                errors_all[label] = errors
        return fvecs_all, errors_all


def compile_featurizer(features_dict):
    """
    From a dict following the feature format specified in the `native` module,
    remove features specifying `compute=custom` and put it into a seperate dict. 
    
    Combine the native and custom featurizers into one dict.
    """
    from .native import compile_native_featurizer
    from .custom_features import compile_custom_featurizer
    native_features = {}
    custom_features_dict = {}
    if (features := features_dict.get("features")) is None:
        raise ValueError("Expected `features` key in the provided dict.")
    for featname, feature_params in list(features.items()):
        if feature_params.get("compute") == "custom":
            custom_features_dict[featname] = feature_params
        else:
            native_features[featname] = feature_params
    native, errors_native = compile_native_featurizer({
        "features": native_features,
        "residue_groups": features_dict.get("residue_groups"),
        "motif_frequencies": features_dict.get("motif_frequencies"),
        "aa_frequencies": features_dict.get("aa_frequencies")
    })
    custom, errors_custom = compile_custom_featurizer(custom_features_dict)
    return_value = native
    return_value.update(custom)
    return_errors = errors_native
    return_errors.update(errors_custom)
    return return_value, return_errors