# SplitClusterTest.jl 


[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://szcf-weiya.github.io/SplitClusterTest.jl/dev) [![codecov](https://codecov.io/gh/szcf-weiya/SplitClusterTest.jl/graph/badge.svg?token=dsRMZFM1q5)](https://codecov.io/gh/szcf-weiya/SplitClusterTest.jl) [![CI](https://github.com/szcf-weiya/SplitClusterTest.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/szcf-weiya/SplitClusterTest.jl/actions/workflows/ci.yml)

Julia package for

> Wang, L., Lin, Y., & Zhao, H. (2024). False Discovery Rate Control via Data Splitting for Testing-after-Clustering (arXiv:2410.06451). arXiv. <https://doi.org/10.48550/arXiv.2410.06451>
>

The proposed approach addresses the double-dipping issue in testing-after-clustering tasks, particularly in single-cell data analysis, where the same data is used both for clustering (to identify cell types) and for testing (to select differentially expressed genes), which can inflate false positives.

![dd](https://github.com/user-attachments/assets/e5383503-2e4d-45d0-adff-77f3a0f82899)

## :arrow_right: See also

- R package: <https://github.com/szcf-weiya/SplitClusterTest>
- For the comparison between data splitting and data fission, check <https://github.com/szcf-weiya/fission_vs_splitting>.
