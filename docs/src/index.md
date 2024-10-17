# SplitClusterTest.jl Documentation

> Wang, L., Lin, Y., & Zhao, H. (2024). False Discovery Rate Control via Data Splitting for Testing-after-Clustering (arXiv:2410.06451). arXiv. <https://doi.org/10.48550/arXiv.2410.06451>
>


Testing for differences in features between clusters in various applications often leads to inflated false positives when practitioners use the same dataset to identify clusters and then test features, an issue commonly known as “double dipping”. 

To address this challenge, inspired by data-splitting strategies for controlling the false discovery rate (FDR) in regressions ([Dai et al., 2023](https://www.tandfonline.com/doi/abs/10.1080/01621459.2022.2060113)), we present a novel method that applies data-splitting to control FDR while maintaining high power in unsupervised clustering. 

We first divide the dataset into two halves, then apply the conventional testing-after-clustering procedure to each half separately and combine the resulting test statistics to form a new statistic for each feature. The new statistic can help control the FDR due to its property of having a sampling distribution that is symmetric around zero for any null feature.
