using Test
using SplitClusterTest


n = 1000
p = 2000
prop = 0.1
truth = 1:Int(p * prop)

@testset "iid Gaussian case" begin
    x, cl = gen_data_normal(n, p, 1.0; prop_imp = prop)
    res = ds(x, q = 0.05, plt = false, ret_ms = false)    
    fdr, f1, prec, recall = calc_acc(res, truth)
    @test 0 <= fdr < 0.2 # although it is 0.05, be slightly loose 
    @test 0.5 <= f1 <= 1.0

    res = mds(x, M = 10, q = 0.05)
    fdr, f1, prec, recall = calc_acc(res, truth)
    @test 0 <= fdr < 0.2
    @test 0.5 <= f1 <= 1.0
end

@testset "correlated Gaussian case" begin
    x, _ = gen_data_normal(n, p, 1.0; prop_imp = prop, ρ = 0.5, corr_structure = "ar1")
    res = ds(x, q = 0.05, plt = false, ret_ms = false)
    fdr, f1, prec, recall = calc_acc(res, truth)
    @test 0 <= fdr < 0.2
    @test 0.5 <= f1 <= 1.0

    res = mds(x, M = 10, q = 0.05)
    fdr, f1, prec, recall = calc_acc(res, truth)
    @test 0 <= fdr < 0.2
    @test 0.5 <= f1 <= 1.0
end

@testset "iid Poisson case" begin
    x = gen_data_pois(n = n, p = p, prop_imp = 0.1, ρ = 0.0, δ = 1.0)
    res = ds(x, q = 0.05, plt = false, ret_ms = false)
    fdr, f1, prec, recall = calc_acc(res, truth)
    @test 0 <= fdr < 0.2
    @test 0.5 <= f1 <= 1.0

    res = mds(x, M = 10, q = 0.05)
    fdr, f1, prec, recall = calc_acc(res, truth)
    @test 0 <= fdr < 0.2
    @test 0.5 <= f1 <= 1.0

    x = gen_data_pois(n = n, p = p, prop_imp = 0.1, ρ = 0.0, δ = 1.0, sigma = 0.1)
    res = ds(x, q = 0.05, plt = false, ret_ms = false)
    fdr, f1, prec, recall = calc_acc(res, truth)
    @test 0 <= fdr < 0.2
    @test 0.5 <= f1 <= 1.0
end
