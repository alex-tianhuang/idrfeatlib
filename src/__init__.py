"""
Another feature analysis library.

Hopefully I finally learned how to write a usable one this time.

Jan 6, 2025
Tian Hao Huang (@tianh) 
"""
import typing

class Featurizer:
    """
    An object that holds a dict of (sequence -> float) functions, i.e. `features`.
    
    Call `Featurizer.featurize` to compute the feature vector of a single sequence.
    """
    _funcs: typing.Dict[str, typing.Callable[..., float]]
    def __init__(self, funcs: typing.Dict[str, typing.Callable[..., float]]):
        """`funcs` should be a dict of (sequence -> float) functions."""
        self._funcs = funcs
    def featurize(self, sequence: str) -> "FeatureVector":
        """Compute the feature vector of a single sequence."""
        return_value = {}
        for featname, func in self._funcs.items():
            return_value[featname] = func(sequence)
        return FeatureVector(return_value)
    def add_feature_function(self, featname: str, func: typing.Callable[[str], float]):
        """Add a feature function."""
        self._funcs[featname] = func
    
class FeatureVector:
    """A dict with arithmetic ops. The `as_dict` attribute contains the underlying data."""
    as_dict: typing.Dict[str, float]
    def __init__(self, data: typing.Dict[str, float]):
        self.as_dict = data
    def __add__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError("cannot add %s to FeatureVector" % type(other).__name__)
        return_value = {}
        for featname, value in self.as_dict.items():
            return_value[featname] = value + other.as_dict[featname]
        return FeatureVector(return_value)
    def __iadd__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError("cannot add %s to FeatureVector" % type(other).__name__)
        for featname in self.as_dict.keys():
            self.as_dict[featname] += other.as_dict[featname]
    def __sub__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError("cannot subtract %s from FeatureVector" % type(other).__name__)
        return_value = {}
        for featname, value in self.as_dict.items():
            return_value[featname] = value - other.as_dict[featname]
        return FeatureVector(return_value)
    def __isub__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError("cannot subtract %s from FeatureVector" % type(other).__name__)
        for featname in self.as_dict.keys():
            self.as_dict[featname] -= other.as_dict[featname]
    def __mul__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError("cannot multiply FeatureVector with %s" % type(other).__name__)
        return_value = {}
        for featname, value in self.as_dict.items():
            return_value[featname] = value * other.as_dict[featname]
        return FeatureVector(return_value)
    def __imul__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError("cannot multiply FeatureVector with %s" % type(other).__name__)
        for featname in self.as_dict.keys():
            self.as_dict[featname] *= other.as_dict[featname]
    def __truediv__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError("cannot divide FeatureVector with %s" % type(other).__name__)
        return_value = {}
        for featname, value in self.as_dict.items():
            return_value[featname] = value / other.as_dict[featname]
        return FeatureVector(return_value)
    def __itruediv__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError("cannot divide FeatureVector with %s" % type(other).__name__)
        for featname in self.as_dict.keys():
            self.as_dict[featname] /= other.as_dict[featname]
    def map_values(self, func):
        """Applies `func` to values of the feature vector."""
        return_value = {
            featname: func(value) for featname, value in self.as_dict.items()
        }
        return FeatureVector(return_value)
    @staticmethod
    def cmv(feature_vectors: typing.Iterable["FeatureVector"]):
        """
        Returns a Count, Mean, Variance (hence the name `cmv`) tuple from the iterable of feature vectors.

        Variance is population variance (i.e. uses `N` instead of `N-1` in the denominator.)
        """
        count = {}
        sum = {}
        sum_sqr = {}
        for fvec in feature_vectors:
            for featname, value in fvec.as_dict.items():
                count[featname] = count.get(featname, 0) + 1
                sum[featname] = sum.get(featname, 0) + value
                sum_sqr[featname] = sum_sqr.get(featname, 0) + value * value
        mean = FeatureVector(sum) / FeatureVector(count)
        var = FeatureVector(sum_sqr) / FeatureVector(count) - mean * mean
        return FeatureVector(count), mean, var
