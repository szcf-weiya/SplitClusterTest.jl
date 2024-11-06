var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [SplitClusterTest]\nOrder = [:type, :function]","category":"page"},{"location":"api/#SplitClusterTest.calc_acc-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{T}}} where T<:Real","page":"API","title":"SplitClusterTest.calc_acc","text":"calc_acc(pred::AbstractVector, truth::AbstractVector)\n\nCalculate the accuracy (FDR, F1 score, Precision, Power) for the predicted selection pred given the truth.\n\n\n\n\n\n","category":"method"},{"location":"api/#SplitClusterTest.calc_τ","page":"API","title":"SplitClusterTest.calc_τ","text":"calc_τ(ms::AbstractVector, q::Float64 = 0.05, offset::Int = 1)\n\nCalculate the cutoff of the mirror statistics ms given the nominal FDR level q. It is recommended to take offset = 1 in the numerator, as discussed in the knockoff paper.\n\n\n\n\n\n","category":"function"},{"location":"api/#SplitClusterTest.ds-Tuple{AbstractMatrix}","page":"API","title":"SplitClusterTest.ds","text":"ds(::AbstractMatrix; ...)\n\nSelect with the nominal FDR level q via a single data splitting on data matrix x.\n\nq: nominal FDR level\nsignal_measure: the signal measurement\nret_ms: if true, then return the mirror statistics; otherwise, return the selection set\ntype: if discrete, perform the testing after clustering into two groups;        otherwise, perform the testing along pseduotime after estimating the pseduotime\ncl_method: the function for clustering into two groups (only used if type == discrete)\nti_method: the function for estimating the pseduotime (only used if type != discrete)\noracle_label: if provided (it is nothing by default), the accuracy of the clustering will be calculated. \nkmeans_whiten: whether to perform whitening\nΣ: used for kmeans whitening\n\n\n\n\n\n","category":"method"},{"location":"api/#SplitClusterTest.first_two_PCs-Tuple{AbstractMatrix}","page":"API","title":"SplitClusterTest.first_two_PCs","text":"first_two_PCs(x::AbstractMatrix)\n\nCalculate the first two principal components of data matrix x.\n\n\n\n\n\n","category":"method"},{"location":"api/#SplitClusterTest.gen_data_normal-Tuple{Int64, Int64, Float64}","page":"API","title":"SplitClusterTest.gen_data_normal","text":"gen_data_normal(n::Int, p::Int, δ::Float64; prop_imp = 0.1, corr_structure = \"ind\", ρ = 0.9, sigma = 0, ...)\n\nGenerate n samples with p features from two Gaussian distributions. \n\nprop_imp: the proportion of relevant features\ncorr_structure: the correlation structure, possible choices:\nind: independent\nar1: AR(1) structure\nfixcorr: fixed correlation\nfixcorr_s1_ind: fixed correlation among relevant features, and no correlation between the null features and relevant features\nfixcorr_s1: fixed correlation among relevant features, and maximum correlation between the null features and relevant features such that the correlation matrix is positive definite.\nρ: the correlation coefficient\nsigma: the noise level on the signal strength\n\n\n\n\n\n","category":"method"},{"location":"api/#SplitClusterTest.gen_data_pois-Tuple{AbstractMatrix}","page":"API","title":"SplitClusterTest.gen_data_pois","text":"gen_data_pois(Λ::AbstractMatrix; ρ = 0.5, block_size = 10)\n\nGenerate Poisson samples of the same size Λ, which is the mean values of each element.  The features can be correlated via the Gaussian Copula, whose correlation matrix is block-diagonal AR(1) structure.\n\nρ: the correlation coefficient\nblock_size: the size of each block in the correlation matrix\n\n\n\n\n\n","category":"method"},{"location":"api/#SplitClusterTest.gen_data_pois-Tuple{Int64, Int64, Float64}","page":"API","title":"SplitClusterTest.gen_data_pois","text":"gen_data_pois(n::Int, p::Int, δ::Float64; prop_imp = 0.1, ρ = 0.5, block_size = 10, type = \"discrete\", sigma = 0)\n\nGenerate n samples with p features from Poisson distributions\n\ntype: if discrete, then generate from two Poisson distributions; otherwise, each sample has a different mean vector, which forms the linear pseduotime data.\nprop_imp: the proportion of relevant features\nsigma: the noise level on the signal strength\nρ: the correlation coefficient. If non-zero, take the Gaussian Copula with AR(1) correlation structure with ρ\nblock_size: the block size when construct the correlation matrix, since the Copula of high dimension is computational expensive.\n\n\n\n\n\n","category":"method"},{"location":"api/#SplitClusterTest.mds-Tuple{AbstractMatrix}","page":"API","title":"SplitClusterTest.mds","text":"mds(x::AbstractMatrix; M = 10, ...)\n\nSelect with the nominal FDR level q via M times data splitting on data matrix x. All paramaters except M are passed to ds.\n\n\n\n\n\n","category":"method"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"EditURL = \"../../../examples/two-gaussians.jl\"","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"This section demonstrates the data splitting procedure for selecting relevant features when there exists (or no) cluster structure under the Gaussian setting.","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"using SplitClusterTest\nusing Plots","category":"page"},{"location":"examples/two-gaussians/#Without-cluster-structure","page":"Two Gaussians","title":"Without cluster structure","text":"","category":"section"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"x, cl = gen_data_normal(1000, 2000, 0.0, prop_imp=0.1);\nnothing #hide","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"Plot the first two PCs of X","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"pc1, pc2 = first_two_PCs(x)\nscatter(pc1[cl .== 1], pc2[cl .== 1])\nscatter!(pc1[cl .== 2], pc2[cl .== 2])","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"perform the data splitting procedure for selecting relevant features","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"ms = ds(x, ret_ms = true);\nnothing #hide","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"the mirror statistics are symmetric about zero since all features are null features.","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"histogram(ms, label = \"\")","category":"page"},{"location":"examples/two-gaussians/#With-cluster-structure","page":"Two Gaussians","title":"With cluster structure","text":"","category":"section"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"x, cl = gen_data_normal(1000, 2000, 0.5, prop_imp=0.1);\nnothing #hide","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"Plot the first two PCs of X","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"pc1, pc2 = first_two_PCs(x)\nscatter(pc1[cl .== 1], pc2[cl .== 1])\nscatter!(pc1[cl .== 2], pc2[cl .== 2])","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"the mirror statistics of relevant features tend to be larger and away from null features, where the null features still exhibit a symmetric distribution about zero.","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"ms = ds(x, ret_ms = true);\nnothing #hide","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"Then we can properly take the cutoff to control the FDR, as shown by the red vertical line.","category":"page"},{"location":"examples/two-gaussians/","page":"Two Gaussians","title":"Two Gaussians","text":"τ = calc_τ(ms)\nhistogram(ms, label = \"\")\nPlots.vline!([τ], label = \"\", lw = 3)","category":"page"},{"location":"examples/cont-poissons/","page":"Continuous Poissons (Linear Pseduotime)","title":"Continuous Poissons (Linear Pseduotime)","text":"EditURL = \"../../../examples/cont-poissons.jl\"","category":"page"},{"location":"examples/cont-poissons/","page":"Continuous Poissons (Linear Pseduotime)","title":"Continuous Poissons (Linear Pseduotime)","text":"This section demonstrates the data splitting procedure for selecting relevant features when there exists latent linear pseduotime under the Poisson setting.","category":"page"},{"location":"examples/cont-poissons/","page":"Continuous Poissons (Linear Pseduotime)","title":"Continuous Poissons (Linear Pseduotime)","text":"using SplitClusterTest\nusing Plots\n\n\nx, cl = gen_data_pois(1000, 2000, 0.5, prop_imp=0.1, type = \"continuous\")","category":"page"},{"location":"examples/cont-poissons/","page":"Continuous Poissons (Linear Pseduotime)","title":"Continuous Poissons (Linear Pseduotime)","text":"Plot the first two PCs of X, and color each point by the pseduotime variable cl","category":"page"},{"location":"examples/cont-poissons/","page":"Continuous Poissons (Linear Pseduotime)","title":"Continuous Poissons (Linear Pseduotime)","text":"pc1, pc2 = first_two_PCs(x)\nscatter(pc1, pc2, marker_z = cl, label = \"\")","category":"page"},{"location":"examples/cont-poissons/","page":"Continuous Poissons (Linear Pseduotime)","title":"Continuous Poissons (Linear Pseduotime)","text":"Adopt the data splitting procedure to select the relevant features.","category":"page"},{"location":"examples/cont-poissons/","page":"Continuous Poissons (Linear Pseduotime)","title":"Continuous Poissons (Linear Pseduotime)","text":"ms = ds(x, ret_ms = true, type = \"continuous\");\nτ = calc_τ(ms)","category":"page"},{"location":"examples/cont-poissons/","page":"Continuous Poissons (Linear Pseduotime)","title":"Continuous Poissons (Linear Pseduotime)","text":"the mirror statistics of relevant features tend to be larger and away from null features, where the null features still exhibit a symmetric distribution about zero. Then we can properly take the cutoff to control the FDR, as shown by the red vertical line.","category":"page"},{"location":"examples/cont-poissons/","page":"Continuous Poissons (Linear Pseduotime)","title":"Continuous Poissons (Linear Pseduotime)","text":"histogram(ms, label = \"\")\nPlots.vline!([τ], label = \"\", lw = 3)","category":"page"},{"location":"#SplitClusterTest.jl-Documentation","page":"Home","title":"SplitClusterTest.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Wang, L., Lin, Y., & Zhao, H. (2024). False Discovery Rate Control via Data Splitting for Testing-after-Clustering (arXiv:2410.06451). arXiv. https://doi.org/10.48550/arXiv.2410.06451","category":"page"},{"location":"","page":"Home","title":"Home","text":"Testing for differences in features between clusters in various applications often leads to inflated false positives when practitioners use the same dataset to identify clusters and then test features, an issue commonly known as “double dipping”. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: dd)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: xkcd)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The xkcd-style cartoon is drawn with the help of R package xkcd","category":"page"},{"location":"","page":"Home","title":"Home","text":"To address this challenge, inspired by data-splitting strategies for controlling the false discovery rate (FDR) in regressions (Dai et al., 2023), we present a novel method that applies data-splitting to control FDR while maintaining high power in unsupervised clustering. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"We first divide the dataset into two halves, then apply the conventional testing-after-clustering procedure to each half separately and combine the resulting test statistics to form a new statistic for each feature. The new statistic can help control the FDR due to its property of having a sampling distribution that is symmetric around zero for any null feature.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: mds)","category":"page"},{"location":"examples/two-poissons/","page":"Two Poissons","title":"Two Poissons","text":"EditURL = \"../../../examples/two-poissons.jl\"","category":"page"},{"location":"examples/two-poissons/","page":"Two Poissons","title":"Two Poissons","text":"This section demonstrates the data splitting procedure for selecting relevant features when there exists cluster structure under the Poisson setting.","category":"page"},{"location":"examples/two-poissons/","page":"Two Poissons","title":"Two Poissons","text":"using SplitClusterTest\nusing Plots\n\n\nx, cl = gen_data_pois(1000, 2000, 0.5, prop_imp=0.1, type = \"discrete\")","category":"page"},{"location":"examples/two-poissons/","page":"Two Poissons","title":"Two Poissons","text":"Plot the first two PCs of X","category":"page"},{"location":"examples/two-poissons/","page":"Two Poissons","title":"Two Poissons","text":"pc1, pc2 = first_two_PCs(x)\nscatter(pc1[cl .== 0], pc2[cl .== 0])\nscatter!(pc1[cl .== 1], pc2[cl .== 1])","category":"page"},{"location":"examples/two-poissons/","page":"Two Poissons","title":"Two Poissons","text":"Adopt the data splitting procedure to select the relevant features.","category":"page"},{"location":"examples/two-poissons/","page":"Two Poissons","title":"Two Poissons","text":"ms = ds(x, ret_ms = true);\nτ = calc_τ(ms)","category":"page"},{"location":"examples/two-poissons/","page":"Two Poissons","title":"Two Poissons","text":"the mirror statistics of relevant features tend to be larger and away from null features, where the null features still exhibit a symmetric distribution about zero. Then we can properly take the cutoff to control the FDR, as shown by the red vertical line.","category":"page"},{"location":"examples/two-poissons/","page":"Two Poissons","title":"Two Poissons","text":"histogram(ms, label = \"\")\nPlots.vline!([τ], label = \"\", lw = 3)","category":"page"}]
}
