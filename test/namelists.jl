using MCTomoTools
using Test

all_lowercase(s) = s == lowercase(s)

@testset "NameLists" begin
    @testset "read_namelist" begin
        file = joinpath(@__DIR__, "data", "sample_MCTomo_input.inp")

        @testset "No types" begin
            nml = MCTomoTools.NameLists.read_namelist(file)
            @test all(key -> haskey(nml, key),
                ("grid_settings", "mcmc_settings", "likelihood_settings"))
            @test nml["grid_settings"]["grid"]["xmin"] == 460
            @test nml["grid_settings"]["grid"]["xmax"] == 473
            @test nml["grid_settings"]["grid"]["nx"] == 101
            @test nml["likelihood_settings"]["like_set"]["datatype"] == 1
            @test nml["likelihood_settings"]["like_set"]["bsources_file"] == "bsources.dat"
            @test nml["likelihood_settings"]["like_set"]["dphasevel"] == 1e-3
            @test nml["likelihood_settings"]["like_set"]["tol"] == 1e-6
            @test nml["mcmc_settings"]["mcmc_set"]["nsamples"] == 1_500_000
            @test nml["mcmc_settings"]["mcmc_set"]["resume"] == 0
            @test nml["mcmc_settings"]["mcmc_set"]["bn1_max"] == 0.55
        end

        @testset "With types" begin
            types = Dict("mcmc_settings" => Dict("mcmc_set" => Dict("pd" => Float32)))
            nml = MCTomoTools.NameLists.read_namelist(file; types=types)
            @test nml["mcmc_settings"]["mcmc_set"]["pd"] isa Float32
        end

        @testset "Lowercase names" begin
            nml = MCTomoTools.NameLists.read_namelist(file)
            for (key, val) in nml
                @test all_lowercase(key)
                for (key, val) in val
                    @test all_lowercase(key)
                    for key in keys(val)
                        @test all_lowercase(key)
                    end
                end
            end
        end
    end
end
