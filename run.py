# -*-coding:utf-8-*-
# Author: Kang Liu <kangliu@stu.xidian.edu.cn>

from preprocessing import preprocessing
from elcnv import ELCNV

from sklearn.cluster import KMeans
from scipy import stats
import pandas as pd
import numpy as np
import warnings


def calculating_CN(mode, CNVRD, CNVtype):
    CN = np.full(len(CNVtype), 0)
    index = CNVtype == 1
    lossRD = CNVRD[index]
    if len(lossRD) > 2:
        data = np.c_[lossRD, lossRD]
        del_type = KMeans(n_clusters=2, random_state=9).fit_predict(data)
        CNVtype[index] = del_type
        if np.mean(lossRD[del_type == 0]) < np.mean(lossRD[del_type == 1]):
            homoRD = np.mean(lossRD[del_type == 0])
            hemiRD = np.mean(lossRD[del_type == 1])
            for i in range(len(CN)):
                if CNVtype[i] == 0:
                    CN[i] = 0
                elif CNVtype[i] == 1:
                    CN[i] = 1
        else:
            hemiRD = np.mean(lossRD[del_type == 0])
            homoRD = np.mean(lossRD[del_type == 1])
            for i in range(len(CN)):
                if CNVtype[i] == 1:
                    CN[i] = 0
                elif CNVtype[i] == 0:
                    CN[i] = 1
        purity = 2 * (homoRD - hemiRD) / (homoRD - 2 * hemiRD)

        for i in range(len(CNVtype)):
            if CNVtype[i] == 2:
                CN[i] = int(2 * CNVRD[i] / (mode * purity) - 2 * (1 - purity) / purity)
    return CN


def combiningCNV(seg_chr, seg_start, seg_end, seg_count, labels, mode):
    def _func(x):
        if x == 2:
            return "duplication"
        else:
            return "deletion"

    index = labels == 1
    CNV_chr = seg_chr[index]
    CNVstart = seg_start[index]
    CNVend = seg_end[index]
    CNVRD = seg_count[index, 0]

    type = np.full(len(CNVRD), 1)
    for i in range(len(CNVRD)):
        if CNVRD[i] > mode:  # "duplication"
            type[i] = 2

    for i in range(len(CNVRD) - 1):
        if CNVend[i] + 1 == CNVstart[i + 1] and type[i] == type[i + 1]:
            CNVstart[i + 1] = CNVstart[i]
            type[i] = 0

    index = type != 0
    CNVRD = CNVRD[index]
    CNV_chr = CNV_chr[index]
    CNVstart = CNVstart[index]
    CNVend = CNVend[index]
    CNVtype = type[index]
    CNVtype = [_func(i) for i in CNVtype]
    return CNV_chr, CNVstart, CNVend, CNVRD, CNVtype


def save_result(CNV_chr, CNVstart, CNVend, CNVRD, CNVtype, bam_path, output_dir=None):
    df = pd.DataFrame()
    df['chr'] = CNV_chr
    df["start"] = CNVstart
    df["end"] = CNVend
    df["type"] = CNVtype
    df["RD"] = CNVRD
    file_name = bam_path.split("/")[-1]
    if output_dir is None:
        df.to_csv(file_name + ".txt", index=False)
    else:
        df.to_csv(output_dir + "/" + file_name + ".txt", index=False)


def sta_score_realdata(groudtruth_path, result_start, result_end, result_type):
    ground_truth = pd.read_table(groudtruth_path)
    # truth_type = ground_truth["variant type"].apply(lambda x: _func(x)).tolist()
    truth_type = ground_truth["variant type"].tolist()
    truth_start = ground_truth['start'].tolist()
    truth_end = ground_truth['stop'].tolist()

    count = 0
    for i in range(len(result_type)):
        for j in range(len(truth_type)):
            if truth_start[j] <= result_start[i] <= truth_end[j] and truth_type[j] == result_type[i]:
                if result_end[i] <= truth_end[j]:
                    count += (result_end[i] - result_start[i] + 1)

                elif result_end[i] > truth_end[j]:
                    count += (truth_end[j] - result_start[i] + 1)

            elif truth_start[j] >= result_start[i] and truth_type[j] == result_type[i]:
                if truth_start[j] <= result_end[i] <= truth_end[j]:
                    count += (result_end[i] - truth_start[j] + 1)

                elif result_end[i] >= truth_end[j]:
                    count += (truth_end[j] - truth_start[j] + 1)

    result_count = 0
    for i in range(len(result_start)):
        result_count += (result_end[i] - result_start[i] + 1)

    truth_count = 0
    for i in range(len(truth_start)):
        truth_count += (truth_end[i] - truth_start[i] + 1)

    # print(count, result_count, truth_count)
    if result_count == 0:
        precision = 0
    else:
        precision = count / result_count
    sensitivity = count / truth_count
    print("ans =", precision, sensitivity)

    return [precision, sensitivity]


def main(bam_path, fa_path, bin_size=1000, output_dir=None, gt_path=None, cbs_imp='python', ncol=50,
         scores_comb=None, is_require_X=False, bandwidth=1.0):
    """
    Parameters
    ----------
    bam_path : str
        Local path of the *.bam file (i.e., sequenced sample).

    fa_path : str
        Local path of the *.fasta file or the *.fa file (i.e., reference genome).

    bin_size : int, optional (default=1000)
        The bin size.

    output_dir : str, optional (default=None)
        Local path for saving experimental results. If the output_dir is None, it is directly
        saved in the path where the code is located.

    gt_path : str, optional (default=None)
        Local path of the ground truth of the sequenced sample.

    cbs_imp: str, optional (default='python')
        The implementation of CBS algorithm. In addition to "python", cbs_imp can also be "R".

    ncol : int, optional (default=50)
        The number of  partitions to CBS in R.

    scores_comb : list, optional (default=["lscp"])
        List of methods are used to merge outlier score vectors in ELCNV. Besides, You can choose the
        following strategies (i.e., ["voting", "maximum", "lscp", "averaging"]).

    is_require_X : bool, optional (default=False)
        When calculating the probability density of the outlier score vector for all segments in NPAT,
        whether the original data X is required.

    bandwidth : float, optional (default=1.0)
        When calculating the probability density of the outlier score vector for all segments in NPAT,
        the bandwidth of the kernel.

    Returns
    ----------

    """

    warnings.filterwarnings("ignore")

    # Preprocessing
    all_chr, all_start, all_end, all_rd, mode = preprocessing(bam_path, fa_path, bin_size=bin_size,
                                                              cbs_imp=cbs_imp, ncol=ncol)

    # ELCNV
    elcnv = ELCNV(scores_comb=scores_comb, is_require_X=is_require_X, bandwidth=bandwidth)
    # Fit individual detectors and predict whether a particular segment is an outlier or not.
    # 0 stands for inlier and 1 for outlier(CNV).
    labels = elcnv.fit_predict(all_rd)

    # Statistics of experimental results
    precision, sensitivity = [], []
    for i in range(labels.shape[1]):
        CNV_chr, CNVstart, CNVend, CNVRD, CNVtype = combiningCNV(all_chr, all_start, all_end, all_rd, labels[:, i],
                                                                 mode)

        # save results
        save_result(CNV_chr, CNVstart, CNVend, CNVRD, CNVtype, bam_path, output_dir)

        if gt_path:
            # statistic the performance(i.e., precision and sensitivity)
            temp_ans = sta_score_realdata(gt_path, CNVstart, CNVend, CNVtype)
            precision.append(temp_ans[0])
            sensitivity.append(temp_ans[1])
    precision = np.round(precision, 2)
    sensitivity = np.round(sensitivity, 2)

    # print results
    for p, s in zip(precision, sensitivity):
        print(f"precision={p:.2f}, recall={s:.2f}, f1-score={stats.hmean((p, s)):.2f}")


if __name__ == '__main__':
    # Local path of the *.bam file
    bam_path = r"/home/mk422/Documents/Code/Python/genetic_analysis/dataset/real_data/NA19238.chrom21.SLX.maq.SRP000032.2009_07.bam"

    # Local path of the *.fasta file or the *.fa file
    fa_path = r"data/chr21.fa"

    # Local path of the ground truth of the sequenced sample.
    gt_path = r"data/NA19238.gt"

    # parameter setting of the preprocessing
    bin_size = 1000  # the bin size ('1000' by default)

    main(bam_path, fa_path, gt_path=gt_path, bin_size=bin_size, cbs_imp='python')