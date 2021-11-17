### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 68a5b8b8-6972-4b03-ad9b-aa0e32c7d426
using Pkg; Pkg.activate(".")

# ╔═╡ caef09cc-0e00-11ec-1753-d7e117eb8c20
begin
	using FITSIO
	using DataFrames
	using LazyArrays, StructArrays
	using SpecialPolynomials
end

# ╔═╡ aed1a891-a5bb-469b-9771-43cb0945d214
using Gumbo

# ╔═╡ fecb4b5e-9046-4cf2-8b82-8f34b09dd662
using ThreadsX

# ╔═╡ 10393148-b9b0-44a0-9b47-c6780318316b
using Statistics: mean

# ╔═╡ 81ab1954-a65e-4633-a340-2e707e5ce879
using FillArrays

# ╔═╡ 4d1ad1b6-788d-4607-b4c3-b37a4cd37e3e
using StaticArrays

# ╔═╡ 05df8c3f-f679-489f-9e51-e5fad58d9a55
using Polynomials

# ╔═╡ ae366e94-e87e-4dec-8424-ae0a4d762ef9
begin
	using PlutoUI, PlutoTest, PlutoTeachingTools
	using BenchmarkTools
	using Profile,  ProfileSVG, FlameGraphs
	using Plots # just needed for colorant when color FlameGraphs
	using LinearAlgebra, PDMats # used by periodogram.jl
	using Random
	Random.seed!(123)
	using FLoops
end;

# ╔═╡ 359c4661-eced-4645-af6e-d862850ab341
md"# Importing Functions From Local File "

# ╔═╡ ae332935-fa0e-48c3-986e-54998dd3ce72
begin
	funcs = ingredients("./src/Support Functions.jl")
	import .funcs: fit_lines_v0_serial, fit_lines_v0_parallel, closest_index, gaussian, gauss_hermite_basis, SpectrumModel, AbsorptionLine, gh_polynomials,  BlazeModel, fit_blaze_model, fit_blaze_model_v0, fit_devs, loss_devs
end

# ╔═╡ 39311305-5492-441b-8346-23d783a04220
md"# Read input data"

# ╔═╡ e624fdcd-df68-4cb7-811a-80154f624a47
begin
	data_path = ""  #Data is stored in same directory as this file
	filename = "neidL1_20210305T170410.fits"   # Just one example
	f = FITS(joinpath(data_path,filename))
	order_idx = 40      # Arbitrarily picking just one order for this example
	pix_fit = 1025:8096  # Exclude pixels near edge of detector where there are issues that complicate things and we don't want to be distracted by



	#reading in wavelength, flux, and variance data
	λ_raw = read(f["SCIWAVE"],:,order_idx)
	flux_raw = read(f["SCIFLUX"],:,order_idx)
	var_raw = read(f["SCIVAR"],:,order_idx)


	#cutting out nans from dataset
	mask_nonans = findall(.!(isnan.(λ_raw) .| isnan.(flux_raw) .| isnan.(var_raw) ))
	λ = λ_raw[mask_nonans]
	flux = flux_raw[mask_nonans]
	var = var_raw[mask_nonans]
end

# ╔═╡ fef84837-477c-4281-8794-c2852f4070eb
md"### Making sure data is read in properly, using tests on the data."

# ╔═╡ 7ce749a5-226a-4207-8984-9b2fad7560ff
#making sure the data is read in correctly--all lambdas have a corresponding flux and variance
begin
	@test length(λ) > 0
	@test length(flux) == length(λ) == length(var)
end

# ╔═╡ 1d610ca0-6981-48ff-b9d8-76c45807b5e7
md"""
# Removing Background Diffraction Pattern
"""

# ╔═╡ 9b91c66f-efa9-4bb5-b0fd-3ebf20a50316
md"""
##### When light passes through a diffraction grating, it takes the shape of a Blaze function, so the observations follow a Blaze shape. We remove the Blaze background here by fitting a model (called Blaze 2) and removing that from the data.
"""

# ╔═╡ ec5afb00-27ae-498f-9a2f-b07afee4e71b
begin  # Removing Blaze background by dividing observations by blaze function
	
	
	mask_fit = pix_fit[findall(.!(isnan.(λ_raw[pix_fit]) .| isnan.(flux_raw[pix_fit]) .| isnan.(var_raw[pix_fit]) ))]
	blaze_model0 =  fit_blaze_model(1:length(λ),flux,var,order=8, mask=mask_fit)
	pix_gt_100p_it1 = flux[pix_fit]./blaze_model0.(pix_fit) .>= 1.0
	
	blaze_model1 =  fit_blaze_model( (1:length(λ)), flux, var,  order=8, mask=pix_fit[pix_gt_100p_it1])
	pix_gt_100p_it2 = flux[pix_fit]./blaze_model1.(pix_fit) .>= 1.0
	
	blaze_model2 =  fit_blaze_model( (1:length(λ)), flux, var,  order=8, mask=pix_fit[pix_gt_100p_it2])
	
	
	local plt = plot()
	plot!(plt,λ[pix_fit],flux[pix_fit], label="Observation",color=:grey)
	plot!(λ[pix_fit],blaze_model2.(pix_fit), label="Blaze 2",color=:blue)
	xlabel!("λ (Å)")
	ylabel!("Flux")

	
	
end

# ╔═╡ cf0842fc-f1b1-416f-be46-c9bc4ceae2be
md"""
##### When the Blaze background is removed, the data takes this shape.
"""

# ╔═╡ c28ee2aa-d3ee-4688-b9a2-cb6b04f31de8
begin
	pix_fit1 = closest_index(λ, 4550):closest_index(λ, 4560)
	local plt = plot(legend=:bottomright)
	#plot!(plt,λ[pix_fit],flux[pix_fit]./blaze_model0.(pix_fit), label="Observation/Blaze 0",color=:red)
	#plot!(plt,λ[pix_fit],flux[pix_fit]./blaze_model1.(pix_fit), label="Observation/Blaze 1",color=:green)
	#plot!(plt,λ[pix_fit1],flux[pix_fit1]./blaze_model2.(pix_fit1), label="Observation/Blaze 2", color=:blue)
	plot!(plt,λ[pix_fit],flux[pix_fit]./blaze_model2.(pix_fit), label="Observation/Blaze 2", color=:blue)

	xlabel!("λ (Å)")
	ylabel!("Normalized Flux")
end

# ╔═╡ 9909c5fb-3e04-4ecf-9399-ccbdf665d35c
md"""
##### Picking a small range of wavelengths. The range is λ = 4569.94 angstroms to λ = 4579.8 angstroms. This data is held in a dataframe df.
"""

# ╔═╡ c9a3550e-ab84-44d4-a935-64e80ed51d63
begin
	
	
	
	pix_plt = closest_index(λ, 4550):closest_index(λ, 4610) #finds the indices that most closely correspond to λ = 4569.94 and λ=4579.8 angstroms. These wavelengths are mostly just arbitrarily chosen, where we found lines within this wavelength band from a list of known absorption lines in the Sun. When we move to parallelized code, we will look at all wavelengths, not just this smaller band.

	# Building a dataframe, where the fluxes have the blaze background removed
df = DataFrame( λ=view(λ,pix_plt), 
				flux=view(flux,pix_plt)./blaze_model2.(pix_plt),
				var =view(var,pix_plt)./(blaze_model2.(pix_plt)).^2
				)
	
end

# ╔═╡ 37309807-cd50-4196-a283-473ee937346a
md"# Fit to Real Solar Absorption Lines"

# ╔═╡ 3137a014-873f-48e5-96d5-49375c3f3ef0
md"""
##### Note: we use a loss function to evaluate the fitting of each line. Loss function is defined as sum[abs(predicted-actual flux)] for λ = λ-3*σ:λ-3*σ
"""

# ╔═╡ f7967636-f51e-4fba-afdb-c8de02da2429
num_gh_orders = 4 #here we can set the order of GH polynomial fits for the rest of the code

# ╔═╡ 26c1319c-8f39-42d8-a3b7-588ac90054d6
begin
	#Found these absorption lines with these standard deviations for each line. 
	#λ_lines = [ 4555.3,4559.95, 4563.23989931, 4563.41910949, 4563.76449209, 4564.17151481, 4564.34024273, 4564.69744168, 4564.82724255, 4565.5196282, 4565.66471622, 4566.23229193, 4566.51949452, 4566.87124681, 4567.40957888, 4568.32947197, 4568.60675666, 4568.77955633, 4569.35589751, 4569.61424051, 4570.02161725, 4570.37846574, 4571.09895901, 4571.43706403, 4571.67596248, 4571.97850945, 4572.27601199, 4572.60413885, 4572.86475639, 4573.80826847, 4573.9777376, 4574.22058447, 4574.47278377, 4574.72227641, 4575.10790085, 4575.54216767, 4575.78807278, 4576.33735793, 4577.17797639, 4577.48371817, 4577.69571021, 4578.03225993, 4578.32482976, 4578.55743589, 4579.05844195, 4579.32999314, 4579.51175328, 4579.67461995, 4579.81960761, 4580.05636844, 4580.41671944, 4580.5881527, 4581.04111864, 4581.20192551, 4582.30743629, 4582.49483022, 4582.83395544, 4583.12727938, 4583.41373121, 4583.83832133, 4584.28107653, 4584.81731555, 4585.08049157, 4585.34095137, 4585.8742966, 4586.22544989, 4586.3712534, 4587.13181633, 4587.72177531, 4588.20292127, 4588.68796051, 4589.01082716, 4589.29840078, 4589.95113751, 4590.7903537, 4591.110535, 4591.39612609, 4592.05460424, 4592.36460681, 4592.65733244, 4593.17249527, 4593.52895394, 4593.92276107, 4594.11902026, 4594.41880008, 4594.63251674, 4594.89747512, 4595.36207176, 4595.59618168, 4596.0598587, 4596.41240259, 4596.57597596, 4596.90736432, 4597.24435651, 4597.38289242, 4597.75193956, 4597.86998627, 4598.12294754, 4598.37433123, 4598.74391477, 4598.99676453, 4599.22776442, 4599.84019949] #length 103

λ_lines_eyeball = [4550.0, 4552.042, 4553.75, 4555.3,4556.76, 4557.2, 4557.4, 4561.341, 4562.63, 4565, 4566.87124681,  4567.69, 4568.3, 4569.5, 4570, 4570.8, 4572.4, 4573.4, 4575.39, 4575.9, 4577.5, 4578.58, 4580, 4581.4, 4581.8, 4582.7, 4584.01, 4585, 4586, 4587.1, 4587.5, 4588.53, 4589.5, 4591.3, 4592.6, 4593.475, 4594, 4595.28, 4596.5, 4597.4, 4599, 4599.5, 4601, 4601.5, 4602, 4603.18, 4604.1, 4605.9, 4606.2, 4606.8, 4607.63, 4608.5, 4609] # has length 53

	λ_lines = λ_lines_eyeball
end;

# ╔═╡ a57b6bcd-9a8b-4f0a-aa98-8a133117785d
begin
	fitted0 = fit_lines_v0_serial(λ_lines,df.λ,df.flux,df.var, order = num_gh_orders)
	fitted_lines = fitted0[1].lines
end


# ╔═╡ 9a06230b-2c9d-43e9-ab7f-00d9b62c44d0
begin
	
	
	losses0 = fitted0[2]
	
end

# ╔═╡ 53975fd7-27e3-4b13-998c-24721dadd0bf
@test length(losses0) == length(fitted_lines) == length(λ_lines) #testing to make sure all the lines have been fitted to, and there is a loss for each fit

# ╔═╡ 0af8099e-0efe-44d2-83e4-abe0d84a900a
md"""
### Regression test for fit
"""

# ╔═╡ 8c3b9f88-4a8b-4d4e-b451-dbf8e2cba69a
@test maximum(losses0) < 15

# ╔═╡ 094e5734-ea75-45d7-b09d-0fe6c6e4c4b5
md"""
### Plotting all fitted lines
"""

# ╔═╡ 27b02672-1c6c-449d-ab7f-5b9dd18bd8a1
begin
	
	fit_windows = fill(0.0,2*length(λ_lines))
	loss_windows = fill(0.0,2*length(λ_lines))
	
	N = length(λ_lines)
	for i in 1:N
		σ_h = 0.04
		fit_windows[2*(i)-1] = λ_lines[i]-fit_devs*σ_h
		loss_windows[2*(i)-1] = λ_lines[i]-loss_devs*σ_h
		
		fit_windows[2*(i)] = λ_lines[i]+fit_devs*σ_h
		loss_windows[2*(i)] = λ_lines[i]+loss_devs*σ_h
	end
	
	
	local plt = plot(legend=:bottomright)
	
	plot!(plt,df.λ,df.flux, label="Observation/Fit")
	plot!(plt,df.λ,fitted0[1](df.λ),label="Model")
	#vline!(λ_lines, label="Line Centers")
	#vline!(fit_windows, label="Fit window")
	#vline!(loss_windows, label="Loss window")
	
	
	xlabel!("λ (Å)")
	ylabel!("Normalized Flux")
end

# ╔═╡ 833904d3-02c1-44c1-b890-219729e9045c
md"""
### Printing out Losses
"""

# ╔═╡ 7913a2aa-6829-4de1-933e-2cb083701b6d
losses0

# ╔═╡ 6208b5ee-eca0-4905-93a9-182992bd3647
fitted_lines

# ╔═╡ 4071e48b-7911-4d90-94d3-746cb9b667ef
md"""
### Printing out Gauss-Hermite coefficients
"""

# ╔═╡ 164aad2d-f4aa-4b29-8fbb-e5ab6758492b
#reading Gauss-Hermite Coefficients from the fitted lines
begin
	
	
	num_lines_found = length(fitted_lines)
	num_gh_coeff = length(fitted_lines[1].gh_coeff)
	gh_s = reshape(zeros(num_lines_found*num_gh_coeff), num_lines_found, num_gh_coeff)
	#gh_s = Array{Float64}(undef, length(lines_found), length(lines_found[1].gh_coeff))
	for i in 1:length(fitted_lines)
		gh = fitted_lines[i].gh_coeff
		for j in 1:length(gh)
			
			gh_s[i,j] = gh[j]
		end
	end
end

# ╔═╡ c4d4523b-e73d-479f-a169-6e40b1f09683
md"""
##### Running tests on the found Gauss-Hermite Coefficients
"""

# ╔═╡ 59aad668-6622-4af9-86c6-0e8f36ef03f7
@test num_lines_found == length(λ_lines) #make sure all lines have been fitted to

# ╔═╡ fe77f6f5-6bd7-47d4-aa49-fa85b555f00e
@test num_gh_coeff == num_gh_orders #make sure all orders exist

# ╔═╡ 71b6349d-f212-4cad-85b8-7d3847dc39cb
#checking if all lines have exactly the same number of gh_coefficients
begin
	all_same_size = 1 #1 as in true, if false, it will be 0
	for i in 1:length(fitted_lines)
		if (length(fitted_lines[i].gh_coeff) != num_gh_coeff)
			all_same_size = 0
		end
	end
	@test all_same_size == 1
end

# ╔═╡ 5252c6b1-d438-416c-a47c-3692be5a2935
with_terminal() do
	for k in 1:Int(length(gh_s)/num_gh_coeff)
		println(gh_s[k, :])	
	end
end

# ╔═╡ c26a6598-3ea3-4f60-b04b-43afc3cbc9bb
md"# Benchmarking"

# ╔═╡ 069f3d91-c287-421a-857d-2fb13fe35fc9
md"# Testing"

# ╔═╡ cc51e5bb-dffb-4a3b-9a21-1c5b426c59c5
begin
	fit_lines_v0_serial(λ_lines,df.λ,df.flux,df.var, order = 
	num_gh_orders)
	fitted_timing_serial = fit_lines_v0_serial(λ_lines,df.λ,df.flux,df.var, order = num_gh_orders)
end

# ╔═╡ 3507a4c3-5602-4c00-9cc1-c82029a69359
begin
	fit_lines_v0_parallel(λ_lines,df.λ,df.flux,df.var, order = num_gh_orders)
	timing_parallel = fit_lines_v0_parallel(λ_lines,df.λ,df.flux,df.var, order = num_gh_orders)
end

# ╔═╡ 061f8147-a6f0-4a1b-9c19-bd65bc37bcee
md"## Packages used"

# ╔═╡ Cell order:
# ╠═68a5b8b8-6972-4b03-ad9b-aa0e32c7d426
# ╠═359c4661-eced-4645-af6e-d862850ab341
# ╠═ae332935-fa0e-48c3-986e-54998dd3ce72
# ╟─39311305-5492-441b-8346-23d783a04220
# ╠═e624fdcd-df68-4cb7-811a-80154f624a47
# ╟─fef84837-477c-4281-8794-c2852f4070eb
# ╠═7ce749a5-226a-4207-8984-9b2fad7560ff
# ╟─1d610ca0-6981-48ff-b9d8-76c45807b5e7
# ╟─9b91c66f-efa9-4bb5-b0fd-3ebf20a50316
# ╠═ec5afb00-27ae-498f-9a2f-b07afee4e71b
# ╟─cf0842fc-f1b1-416f-be46-c9bc4ceae2be
# ╠═c28ee2aa-d3ee-4688-b9a2-cb6b04f31de8
# ╟─9909c5fb-3e04-4ecf-9399-ccbdf665d35c
# ╠═c9a3550e-ab84-44d4-a935-64e80ed51d63
# ╟─37309807-cd50-4196-a283-473ee937346a
# ╟─3137a014-873f-48e5-96d5-49375c3f3ef0
# ╠═f7967636-f51e-4fba-afdb-c8de02da2429
# ╠═26c1319c-8f39-42d8-a3b7-588ac90054d6
# ╠═a57b6bcd-9a8b-4f0a-aa98-8a133117785d
# ╠═9a06230b-2c9d-43e9-ab7f-00d9b62c44d0
# ╠═53975fd7-27e3-4b13-998c-24721dadd0bf
# ╟─0af8099e-0efe-44d2-83e4-abe0d84a900a
# ╠═8c3b9f88-4a8b-4d4e-b451-dbf8e2cba69a
# ╟─094e5734-ea75-45d7-b09d-0fe6c6e4c4b5
# ╠═27b02672-1c6c-449d-ab7f-5b9dd18bd8a1
# ╟─833904d3-02c1-44c1-b890-219729e9045c
# ╠═7913a2aa-6829-4de1-933e-2cb083701b6d
# ╠═6208b5ee-eca0-4905-93a9-182992bd3647
# ╟─4071e48b-7911-4d90-94d3-746cb9b667ef
# ╠═164aad2d-f4aa-4b29-8fbb-e5ab6758492b
# ╟─c4d4523b-e73d-479f-a169-6e40b1f09683
# ╠═59aad668-6622-4af9-86c6-0e8f36ef03f7
# ╠═fe77f6f5-6bd7-47d4-aa49-fa85b555f00e
# ╠═71b6349d-f212-4cad-85b8-7d3847dc39cb
# ╠═5252c6b1-d438-416c-a47c-3692be5a2935
# ╟─c26a6598-3ea3-4f60-b04b-43afc3cbc9bb
# ╟─069f3d91-c287-421a-857d-2fb13fe35fc9
# ╠═cc51e5bb-dffb-4a3b-9a21-1c5b426c59c5
# ╠═3507a4c3-5602-4c00-9cc1-c82029a69359
# ╟─061f8147-a6f0-4a1b-9c19-bd65bc37bcee
# ╠═caef09cc-0e00-11ec-1753-d7e117eb8c20
# ╠═aed1a891-a5bb-469b-9771-43cb0945d214
# ╠═fecb4b5e-9046-4cf2-8b82-8f34b09dd662
# ╠═10393148-b9b0-44a0-9b47-c6780318316b
# ╠═81ab1954-a65e-4633-a340-2e707e5ce879
# ╠═4d1ad1b6-788d-4607-b4c3-b37a4cd37e3e
# ╠═05df8c3f-f679-489f-9e51-e5fad58d9a55
# ╠═ae366e94-e87e-4dec-8424-ae0a4d762ef9
