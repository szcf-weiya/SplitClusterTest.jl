# This section demonstrates the data splitting procedure for selecting relevant features when there exists latent linear pseduotime under the Poisson setting. 

using SplitClusterTest
using Plots


x, cl = gen_data_pois(1000, 2000, 0.5, prop_imp=0.1, type = "continuous")

# Plot the first two PCs of X, and color each point by the pseduotime variable `cl`
pc1, pc2 = first_two_PCs(x)
scatter(pc1, pc2, marker_z = cl, label = "")

# Adopt the data splitting procedure to select the relevant features.
ms = ds(x, ret_ms = true, type = "continuous");
τ = calc_τ(ms)

# the mirror statistics of relevant features tend to be larger and away from null features, 
# where the null features still exhibit a symmetric distribution about zero. 
# Then we can properly take the cutoff to control the FDR, as shown by the red vertical line.
histogram(ms, label = "")
Plots.vline!([τ], label = "", lw = 3)

