# This section demonstrates the data splitting procedure for selecting relevant features when there exists (or no) cluster structure under the Gaussian setting. 

using SplitClusterTest
using Plots


# ## Without cluster structure 
x, cl = gen_data_normal(1000, 2000, 0.0, prop_imp=0.1);

# Plot the first two PCs of X
pc1, pc2 = first_two_PCs(x)
scatter(pc1[cl .== 1], pc2[cl .== 1])
scatter!(pc1[cl .== 2], pc2[cl .== 2])

# perform the data splitting procedure for selecting relevant features
ms = ds(x, ret_ms = true);

# the mirror statistics are symmetric about zero since all features are null features. 
histogram(ms, label = "")


# ## With cluster structure

x, cl = gen_data_normal(1000, 2000, 0.5, prop_imp=0.1);

# Plot the first two PCs of X
pc1, pc2 = first_two_PCs(x)
scatter(pc1[cl .== 1], pc2[cl .== 1])
scatter!(pc1[cl .== 2], pc2[cl .== 2])


# the mirror statistics of relevant features tend to be larger and away from null features, 
# where the null features still exhibit a symmetric distribution about zero. 
ms = ds(x, ret_ms = true);

# Then we can properly take the cutoff to control the FDR, as shown by the red vertical line.
τ = calc_τ(ms)
histogram(ms, label = "")
Plots.vline!([τ], label = "", lw = 3)

