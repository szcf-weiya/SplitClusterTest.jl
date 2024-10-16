# This section demonstrates the distribution of the mirror statistics $\{M_j \}_{j=1}^p$ with or without cluster structure. 

using SplitClusterTest
using Plots

# Without cluster structure, the mirror statistics are symmetric about zero since all features are null features. 
x, cl = gen_data_normal(1000, 2000, 0.5, prop_imp=0.1);
ms = ds(x, method = "tstat", ret_ms = true)
histogram(ms, label = "")


# With cluster structure, the mirror statistics of DE genes tend to be larger and away from null features, 
# where the null features still exhibit a symmetric distribution about zero. 
# Then we can properly take the cutoff to control the FDR, as shown by the red vertical line.
x, cl = gen_data_normal(1000, 2000, 0.0, prop_imp=0.1);
ms = ds(x, method = "tstat", ret_ms = true) # without vertical line
τ = calc_τ(ms)
histogram(ms, label = "")
Plots.vline!([τ], label = "", lw = 3)

