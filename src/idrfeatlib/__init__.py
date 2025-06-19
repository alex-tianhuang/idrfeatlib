"""
An IDR feature analysis library.

Unlike my previous attempts at such a library, this one is very quick and dirty,
not optimized, and only really cares about the "happy path". If you feed
in the wrong type or feed in a strange string as a protein sequence, that's on you
and I will let undefined behaviour take over.

Breakdown
---------
This module contains the `FeatureVector` class, a wrapper around a dict
of floats, with arithmetic operations. To make sure you are always combining
feature vectors with the same shape, use `FeatureVector.shape_eq` to compare
avilable features.

There are three main submodules so far.
1. `featurizer`, which contains the `Featurizer` class.
    This class is a wrapper around a dictionary of function pointers,
    and you can use this class to compute one feature vector or many.
2. `metric`, which contains the `Metric` class.
    This class holds an `origin` and a `weights` vector.
    Norm/distance calculations are done using the formula (x - y) * w
3. `native`, which contains various simple IDR features and also the
    `compile_native_featurizer` function, which returns a dict of function
    pointers, ready for `Featurizer`.

Jan 6, 2025
Tian Hao Huang (@tianh) 
"""
import typing

__all__ = ["FeatureVector"]


class FeatureVector:
    """A dict with arithmetic ops. The `as_dict` attribute contains the underlying data."""

    as_dict: typing.Dict[str, float]

    def __init__(self, data: typing.Dict[str, float]):
        self.as_dict = data

    def __add__(self, other) -> "FeatureVector":
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

    def __sub__(self, other) -> "FeatureVector":
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError(
                "cannot subtract %s from FeatureVector" % type(other).__name__
            )
        return_value = {}
        for featname, value in self.as_dict.items():
            return_value[featname] = value - other.as_dict[featname]
        return FeatureVector(return_value)

    def __isub__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError(
                "cannot subtract %s from FeatureVector" % type(other).__name__
            )
        for featname in self.as_dict.keys():
            self.as_dict[featname] -= other.as_dict[featname]

    def __mul__(self, other) -> "FeatureVector":
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError(
                "cannot multiply FeatureVector with %s" % type(other).__name__
            )
        return_value = {}
        for featname, value in self.as_dict.items():
            return_value[featname] = value * other.as_dict[featname]
        return FeatureVector(return_value)

    def __imul__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError(
                "cannot multiply FeatureVector with %s" % type(other).__name__
            )
        for featname in self.as_dict.keys():
            self.as_dict[featname] *= other.as_dict[featname]

    def __truediv__(self, other) -> "FeatureVector":
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError(
                "cannot divide FeatureVector with %s" % type(other).__name__
            )
        return_value = {}
        for featname, value in self.as_dict.items():
            return_value[featname] = value / other.as_dict[featname]
        return FeatureVector(return_value)

    def __itruediv__(self, other):
        if isinstance(other, dict):
            other = FeatureVector(other)
        if not isinstance(other, FeatureVector):
            raise TypeError(
                "cannot divide FeatureVector with %s" % type(other).__name__
            )
        for featname in self.as_dict.keys():
            self.as_dict[featname] /= other.as_dict[featname]

    def __eq__(self, other):
        if not isinstance(other, FeatureVector):
            raise TypeError(
                "cannot compare FeatureVector with %s" % type(other).__name__
            )
        return self.as_dict == other.as_dict

    def shape_eq(self, other: "FeatureVector"):
        return self.as_dict.keys() == other.as_dict.keys()

    def map_values(self, func) -> "FeatureVector":
        """Applies `func` to values of the feature vector."""
        return_value = {
            featname: func(value) for featname, value in self.as_dict.items()
        }
        return FeatureVector(return_value)

    @staticmethod
    def cmv(
        feature_vectors: typing.Iterable["FeatureVector"],
    ) -> typing.Tuple["FeatureVector", "FeatureVector", "FeatureVector"]:
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

    @staticmethod
    def dump(
        feature_vectors: typing.List[typing.Tuple[typing.Any, "FeatureVector"]],
        path: str,
        label_names: typing.Any = "label",
    ):
        """
        Save labelled feature vectors to csv.

        Usage
        -----
        The "labels" can be strings or ints (i.e. `feature_vectors` can be a list of (key, featvec) pairs),
        and `label_names` is expected to be a string in this case.

        The "labels" can also be n-tuples (i.e. `feature_vectors` can be a list of ((key1, key2, ...), featvec)))
        and `label_names` is expected to be a n-tuple of strings (e.g. ("label1", "label2", etc.))
        """
        import csv

        if not (isinstance(label_names, str) or isinstance(label_names, tuple)):
            raise TypeError("`label_names` must be a string or n-tuple of strings")
        if not feature_vectors:
            if isinstance(label_names, str):
                label_names = (label_names,)
            with open(path, "w") as file:
                writer = csv.DictWriter(file, fieldnames=label_names)
                writer.writeheader()
            return
        featnames = set()
        ordered_featnames = list()
        for _, feature_vector in feature_vectors:
            for featname in feature_vector.as_dict.keys():
                if featname not in featnames:
                    featnames.add(featname)
                    ordered_featnames.append(featname)
        if isinstance(label_names, str):
            colnames = [label_names] + ordered_featnames
        else:
            colnames = list(label_names) + ordered_featnames
        with open(path, "w") as file:
            writer = csv.DictWriter(file, fieldnames=colnames)
            writer.writeheader()
            for label, fvec in feature_vectors:
                if isinstance(label_names, str):
                    assert isinstance(
                        label, str
                    ), "if `label_names` is a string, `feature_vectors` should look like this: [(key, featvec)]"
                    writer.writerow({**fvec.as_dict, label_names: label})
                else:
                    row = {**fvec.as_dict}
                    assert (
                        isinstance(label, tuple) and len(label) == len(label_names)
                    ), "if `label_names` is a tuple, `feature_vectors` should look like this: [((key1, key2, ...), featvec)]"
                    for sublabel, label_name in zip(label, label_names):
                        row[label_name] = sublabel
                    writer.writerow(row)

    @staticmethod
    def load(
        path: str, labels: typing.Any = None
    ) -> typing.Iterator[typing.Tuple[typing.Any, "FeatureVector"]]:
        """
        Load feature vectors from csv (as an iterator).

        Usage
        -----
        The `labels` argument determines what type the iterator will yield.

        In the return type `Iterator[Tuple[X, FeatureVector]]`, I refer to the type `X` as the "left-type".

        If `labels` is None, left-type are strings corresponding to the first column of the csv.
        If `labels` is a string, the left-type are strings corresponding to that column of the csv.
        If `labels` is a tuple of strings, the left-type are string tuples corresponding to those columns of the csv.

        Example
        -------
        fvecs = list(FeatureVector.load("some_fvecs.csv"))
        """
        import csv

        with open(path, "r") as file:
            reader = csv.DictReader(file)
            if labels is None:
                if reader.fieldnames is None:
                    return
                labels = reader.fieldnames[0]
            if isinstance(labels, tuple):
                for row in reader:
                    label_tuple = tuple(row.pop(column) for column in labels)
                    yield (
                        label_tuple,
                        FeatureVector(
                            {
                                featname: float(value)
                                for featname, value in row.items()
                                if value
                            }
                        ),
                    )
            else:
                for row in reader:
                    label = row.pop(labels)
                    yield (
                        label,
                        FeatureVector(
                            {
                                featname: float(value)
                                for featname, value in row.items()
                                if value
                            }
                        ),
                    )

    def __repr__(self) -> str:
        return "FeatureVector(%s)" % repr(self.as_dict)

    def __str__(self) -> str:
        return "FeatureVector(featdim=%d)" % len(self.as_dict)
