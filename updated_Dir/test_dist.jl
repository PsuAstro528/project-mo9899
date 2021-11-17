using Distributed

@everywhere begin
	using Pkg
	#Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true
	Pkg.activate(".")
	#Pkg.instantiate()
	#Pkg.precompile()
	include("./src/Support Functions.jl")
end

@everywhere using .SupportFunctions

@sync @distributed for i in 1:10
	test()
end