# ELCNV
An ensemble learning approach to detect CNVs from NGS data using ELCNV </br>

# Usage
1. API for the entire process
```python
# (1) open the `run.py` and modify the `bam_path` and the `fa_path`
if __name__ == '__main__':
    # Local path of the *.bam file
    bam_path = r"/real_data/NA19238.chrom21.SLX.maq.SRP000032.2009_07.bam"
    
    # Local path of the *.fasta file or the *.fa file
    fa_path = r"data/chr21.fa"
    
    # Local path of the ground truth of the sequenced sample.
    gt_path = r"data/NA19238.gt"
    
    # parameter setting of the preprocessing
    bin_size = 1000  # the bin size ('1000' by default)

    main(bam_path, fa_path, gt_path=gt_path, bin_size=bin_size, cbs_imp='python')
# (2) run the `run.py` 
# (3) If you want to know more detailed introduction, you can check this `main` function in the `run.py`
```

2. API only for the core module of the ELCNV
```python
from elcnv import ELCNV

clf = ELCNV()

# 0 stands for inlier and 1 for outlier(CNV).
# labels with different merging strategies are obtained by the ELCNV and the NPAT method
labels = clf.fit_predict(X) 

# If you only need anomaly score vectors with different merging strategies.
scores = clf.decision_function(X)
```

3. API only for the NPAT method
```python
from npat import NPAT

clf = NPAT()
clf.fit(X)

# 0 stands for inlier and 1 for outlier(CNV).
labels = clf.labels_

```


# Required Dependencies
1. Python 3.8            
    - biopython     1.78
    - combo         0.1.1
    - numpy         1.18.5
    - pandas        1.0.5
    - pysam         0.16.0.1
    - pyod          0.8.4
    - rpy2          3.4.2
    - scikit-learn  0.23.1
    - scipy         1.5.0
2. R 3.4.4
    - DNAcopy
