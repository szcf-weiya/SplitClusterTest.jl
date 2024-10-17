# This section demonstrates the data splitting procedure for selecting relevant features when there exists cluster structure under the Poisson setting. 

using SplitClusterTest
using Plots


x, cl = gen_data_pois(1000, 2000, 0.5, prop_imp=0.1, type = "discrete")

# Plot the first two PCs of X
pc1, pc2 = first_two_PCs(x)
scatter(pc1[cl .== 0], pc2[cl .== 0])
scatter!(pc1[cl .== 1], pc2[cl .== 1])

# Adopt the data splitting procedure to select the relevant features.
ms = ds(x, ret_ms = true);
τ = calc_τ(ms)

# the mirror statistics of relevant features tend to be larger and away from null features, 
# where the null features still exhibit a symmetric distribution about zero. 
# Then we can properly take the cutoff to control the FDR, as shown by the red vertical line.
histogram(ms, label = "")
Plots.vline!([τ], label = "", lw = 3)

