# -*-coding:utf-8-*-
"""A Binary Clustering Model (BCM)
"""
# Author: Kang Liu <kangliu@stu.xidian.edu.cn>

from __future__ import division
from __future__ import print_function

import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.cluster import KMeans
from sklearn.utils.validation import check_array


class BCM(object):
    """
    A new binary clustering model (BCM) without an assumption of
    any distributions.
    This method first employs the kernel density estimation (KDE) method
    with default parameters to estimate the probability density function
    of the outlier score vector, then leverages the 2-means method to
    bi-cluster its probability density. Finally, the cluster with a
    smaller probability density is regarded as CNVs.

    Parameters
    ----------
    X : numpy array of shape (n_samples, n_features), optional (default=None)
        The original data

    bandwidth : float, optional (default=1.0)
        The bandwidth of the kernel.

    is_require_X : bool, optional (default=False)
        When calculating the probability density of the outlier score vector, whether
        the original data X is required.

    Attributes
    ----------
    labels_ : numpy array of shape (n_samples, )
            Binary labels of all merging strategies to indicate whether each segment is a CNV.
            0 stands for inlier and 1 for outlier(CNV).
            
    Examples:
        .. testcode::

            clf = BCM()
            labels = clf.fit_predict(anomaly_scores.reshape(-1, 1))

    """
    def __init__(self, X=None, bandwidth=1.0, is_require_X=False):
        self.X = X
        self.bandwidth = bandwidth
        self.is_require_X = is_require_X
        self.labels_ = None

    def fit(self, anomaly_scores):
        """
        Fit estimator

        Parameters
        ----------
        anomaly_scores : numpy array of shape (n_samples, 1)
            The input samples.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        anomaly_scores = check_array(anomaly_scores)
        kde = KernelDensity(bandwidth=self.bandwidth)

        # The original data X is also used to calculate the probability density
        if self.is_require_X:
            self.X = check_array(self.X)
            anomaly_scores = np.hstack((self.X, anomaly_scores))
        kde.fit(anomaly_scores)
        proba = np.exp(kde.score_samples(anomaly_scores))  # original version
        # note that since the original version added an exponential function, 
        # it will reduce the density difference between inlier and outlier samples.
        # proba = kde.score_samples(anomaly_scores)   # log likelihood

        cluster = KMeans(n_clusters=2, algorithm="elkan")
        _labels = cluster.fit_predict(proba.reshape(-1, 1))
        anomal_label = np.argmin(cluster.cluster_centers_)
        self.labels_ = _labels == anomal_label
        # outlier: 1; inlier: 0
        self.labels_ = self.labels_.astype(int)

    def fit_predict(self, anomaly_scores):
        """
        Fit individual detectors and predict whether a particular segment
        is an outlier or not.

        Parameters
        ----------
        anomaly_scores : numpy array of shape (n_samples, 1)
            The input samples.

        Returns
        -------
        labels : numpy array of shape (n_samples, len(self.scores_comb))
            Binary labels of all merging strategies to indicate whether each segment is a CNV.
            0 stands for inlier and 1 for outlier(CNV).
        """
        if self.labels_ is None:
            self.fit(anomaly_scores)
        return self.labels_
