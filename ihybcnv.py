# -*-coding:utf-8-*-
"""IhybCNV: an Intra-hybrid Approach for CNV Detection from Next-Generation Sequencing Data.
"""
# Author: Kang Liu <kangliu@stu.xidian.edu.cn>

from __future__ import division
from __future__ import print_function

import numpy as np
import statistics
from sklearn.utils import check_array
from sklearn.preprocessing import scale
from pyod.models.lof import LOF
from pyod.models.hbos import HBOS
from pyod.models.cblof import CBLOF
from pyod.models.iforest import IForest
from pyod.models.so_gaal import SO_GAAL
from combo.models.detector_lscp import LSCP

from bcm import BCM


class IhybCNV(object):
    """
    IhybCNV method
    an Intra-hybrid Approach for CNV Detection from Next-Generation Sequencing Data.

    Parameters
    ----------
    scores_comb : list, optional (default=["lscp"])
        List of methods are used to merge outlier score vectors in IhybCNV. Besides, You can choose the
        following strategies (i.e., ["voting", "maximum", "lscp", "averaging"]).

    is_require_X : bool, optional (default=False)
        When calculating the probability density of the outlier score vector for all segments in BCM,
        whether the original data X is required.

    bandwidth : float, optional (default=1.0)
        When calculating the probability density of the outlier score vector for all segments in BCM,
        the bandwidth of the kernel.

    Attributes
    ----------
    labels_ : numpy array of shape (n_samples, len(self.scores_comb))
            Binary labels of all merging strategies to indicate whether each segment is a CNV.
            0 stands for inlier and 1 for outlier(CNV).

    scores_ : numpy array of shape (n_samples, len(self.scores_comb))
            The outlier score vector of all segments using merging strategies.
            Note: when the merging strategies contain voting, its corresponding
            anomaly_score is NaN.

    labels_base_ : numpy array of shape (n_samples, len(detectors))
            Binary labels of individual detectors to indicate whether each segment is a CNV.
            0 stands for inlier and 1 for outlier(CNV).

    scores_base_ : numpy array of shape (n_samples, len(detectors))
            The outlier score vector of individual detectors for all segments.
    """

    def __init__(self, scores_comb=None, is_require_X=False, bandwidth=1.0):
        self.is_require_X = is_require_X
        self.bandwidth = bandwidth
        if scores_comb is None:
            self.scores_comb = ["lscp"]
        else:
            if not isinstance(scores_comb, list):
                raise TypeError("The combination of outlier score vectors need to be a list, "
                                "but what you enter is a %s" % type(scores_comb))
            # Optional schemes of merging outliers scores
            _available_comb = ["voting", "maximum", "lscp", "averaging"]
            for i in range(len(scores_comb)):
                if scores_comb[i].lower() not in _available_comb:
                    raise ValueError("This merging strategy %s is temporarily not supported! "
                                     "You can choose the following strategies." % scores_comb[i], _available_comb)
                else:  # All merging strategies are converted to lowercase
                    scores_comb[i] = scores_comb[i].lower()
            self.scores_comb = scores_comb

        # record results for individual detectors
        self.scores_base_ = None
        self.labels_base_ = None

        # record results for all merging strategies
        self.scores_ = None
        self.labels_ = None

    def fit(self, X):
        """
        Fit individual detectors.

        Parameters
        ----------
        X : numpy array of shape (n_samples, n_features)
            The RD profile of all segments generated after preprocessing.

        Returns
        -------
        self : object
             Fitted estimator.
        """
        X = check_array(X)

        # normalization of all segments with Z-score
        scale_X = scale(X)

        # all base detectors with default parameters
        detectors = [LOF(), SO_GAAL(), IForest(), HBOS(), CBLOF()]

        # record results for individual detectors
        self.scores_base_ = np.zeros((len(scale_X), len(detectors)))
        self.labels_base_ = np.zeros((len(scale_X), len(detectors)))

        # record results for all merging strategies
        self.scores_ = np.zeros((len(scale_X), len(self.scores_comb)))
        self.labels_ = np.zeros((len(scale_X), len(self.scores_comb)))

        for i in range(len(detectors)):
            clf = detectors[i].fit(scale_X)
            self.scores_base_[:, i] = clf.decision_function(scale_X)

            # obtain a series of binary labels using the BCM
            _npat = BCM(X=scale_X, is_require_X=self.is_require_X, bandwidth=self.bandwidth)
            _npat.fit(self.scores_base_[:, i].reshape(-1, 1))
            self.labels_base_[:, i] = _npat.labels_

        # normalization of all outlier score vectors with Z-score
        _scale_score = scale(self.scores_base_)

        for i in range(len(self.scores_comb)):
            if self.scores_comb[i] == "voting":  # majority_vote
                self.scores_[:, i] = np.array([np.nan] * len(scale_X))
                self.labels_[:, i] = np.array([statistics.mode(j) for j in self.labels_base_])

            elif self.scores_comb[i] == "maximum":
                # the maximum of five outlier scores for each segment
                self.scores_[:, i] = np.max(_scale_score, axis=1)

                # obtain binary labels with BCM
                _npat = BCM(X=scale_X, is_require_X=self.is_require_X, bandwidth=self.bandwidth)
                _npat.fit(self.scores_[:, i].reshape(-1, 1))
                self.labels_[:, i] = _npat.labels_

            elif self.scores_comb[i] == "lscp":
                clf = LSCP(detectors, pre_fitted=True)
                clf.fit(scale_X)
                self.scores_[:, i] = clf.decision_function(scale_X)

                # obtain binary labels with the BCM
                _npat = BCM(X=scale_X, is_require_X=self.is_require_X, bandwidth=self.bandwidth)
                _npat.fit(self.scores_[:, i].reshape(-1, 1))
                self.labels_[:, i] = _npat.labels_

            elif self.scores_comb[i] == "averaging":
                self.scores_[:, i] = np.mean(_scale_score, axis=1)

                # obtain binary labels with the BCM
                _npat = BCM(X=scale_X, is_require_X=self.is_require_X, bandwidth=self.bandwidth)
                _npat.fit(self.scores_[:, i].reshape(-1, 1))
                self.labels_[:, i] = _npat.labels_

    def fit_predict(self, X):
        """
        Fit individual detectors and predict whether a particular segment
        is an outlier or not.

        Parameters
        ----------
        X : numpy array of shape (n_samples, n_features)
            The RD profile of all segments generated after preprocessing.

        Returns
        -------
        labels : numpy array of shape (n_samples, len(self.scores_comb))
            Binary labels of all merging strategies to indicate whether each segment is a CNV.
            0 stands for inlier and 1 for outlier(CNV).
        """
        if self.labels_ is None:
            self.fit(X)
        return self.labels_

    def decision_function(self, X):
        """
        Fit individual detectors and return outlier score vectors
        with different merging strategies.

        Parameters
        ----------
        X : numpy array of shape (n_samples, n_features)
            The RD profile of all segments generated after preprocessing.

        Returns
        -------
        scores : numpy array of shape (n_samples, len(self.scores_comb))
            The outlier score vector of all segments using merging strategies.

        Notes
        -----
        when the merging strategies contain voting, its corresponding
        anomaly_score is NaN.
        """
        if self.scores_ is None:
            self.fit(X)
        return self.scores_
