using Distributed

@everywhere begin
	include("./src/Support Functions.jl")
end

@everywhere using .SupportFunctions

@sync @distributed for i in 1:10
	test()
end