from . import FeatureVector
import typing

class Metric:
    """Object holding an `origin` and a `weights` vector."""
    def __init__(self, origin: FeatureVector, weights: FeatureVector) -> None:
        self.origin = origin
        self.weights = weights
    def euclidean_norm_of(self, fvec: FeatureVector) -> float:
        """Compute the euclidean distance from a feature vector to `self.origin`"""
        return self.euclidean_distance_between(fvec, self.origin)
    def euclidean_distance_between(self, fvec_a: FeatureVector, fvec_b: FeatureVector) -> float:
        """Compute the euclidean distance between two feature vectors."""
        z = (fvec_a - fvec_b) * self.weights
        norm = 0
        for coordinate in z.as_dict.values():
            norm += coordinate * coordinate
        return norm
    @staticmethod
    def load(path: str, labels: typing.Any = None, origin_label = "origin", weights_label = "weights"):
        """
        Load an origin and weight feature vector from a csv.

        If no primary id is given, it is assumed to be the first column.

        Example
        -------
        metric = Metric.load("origin-weight.csv")
        """
        fvecs = dict(FeatureVector.load(path, labels=labels))
        origin = fvecs[origin_label]
        weights = fvecs[weights_label]
        return Metric(origin, weights)
    def dump(self, path: str, primary_id = "label", origin_label = "origin", weights_label = "weights"):
        """Write out the origin and weights feature vector to a csv."""
        FeatureVector.dump([
            (origin_label, self.origin),
            (weights_label, self.weights)
        ], path, label_names=primary_id)
        