### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ caef09cc-0e00-11ec-1753-d7e117eb8c20
begin
	using FITSIO
	using DataFrames
	using LazyArrays, StructArrays
	using SpecialPolynomials
end

# ╔═╡ 10393148-b9b0-44a0-9b47-c6780318316b
using Statistics: mean

# ╔═╡ 81ab1954-a65e-4633-a340-2e707e5ce879
using FillArrays

# ╔═╡ 6b4745bf-e13f-4f20-9bda-4dd37662325b
using PDMats

# ╔═╡ 4d1ad1b6-788d-4607-b4c3-b37a4cd37e3e
using StaticArrays

# ╔═╡ 05df8c3f-f679-489f-9e51-e5fad58d9a55
using Polynomials

# ╔═╡ 83cd8a9e-7bdc-4834-b910-8068767ebfac
using PlutoUI

# ╔═╡ e533ddde-2ed7-42b3-a1cd-a1e5738cf3b3
using PlutoTest

# ╔═╡ ce588b8c-4656-4b58-b4e1-9c0ae8b9eefd
using Plots

# ╔═╡ 39311305-5492-441b-8346-23d783a04220
md"## Read input data"

# ╔═╡ e624fdcd-df68-4cb7-811a-80154f624a47
begin
	data_path = ""  # Replace with path to where you store data
	filename = "neidL1_20210305T170410.fits"   # Just one example
end

# ╔═╡ 93b5f3b9-8d7a-43e4-9901-d0f8800b2b28
f = FITS(joinpath(data_path,filename))

# ╔═╡ 79f795c4-a1a6-4e08-a3f1-3bd7bca450ba
order_idx = 40      # Arbitrarily picking just one order for this example

# ╔═╡ 0a132ecb-c759-4448-823a-9c4727ee6916
pix_fit = 1025:8096  # Exclude pixels near edge of detector where there are issues that complicate things and we don't want to be distracted by

# ╔═╡ 1870ae31-79bb-4e5d-a7f1-6202cca94867
begin    # Read in data and exclude any pixels that have issues
	λ = read(f["SCIWAVE"],:,order_idx)
	flux = read(f["SCIFLUX"],:,order_idx)
	var = read(f["SCIVAR"],:,order_idx)
	mask_nonans = findall(.!(isnan.(λ) .| isnan.(flux) .| isnan.(var) ))
	mask_fit = pix_fit[findall(.!(isnan.(λ[pix_fit]) .| isnan.(flux[pix_fit]) .| isnan.(var[pix_fit]) ))]
end

# ╔═╡ 1d610ca0-6981-48ff-b9d8-76c45807b5e7
md"""
## Fit a model to blaze function
"""

# ╔═╡ c9a3550e-ab84-44d4-a935-64e80ed51d63
pix_plt = 3500:4500   # Pick some small range of wavelengths for testing

# ╔═╡ 8d935f79-f2da-4d41-8dad-85cd08197d17
md"## Fit one line (for testing purposes)"

# ╔═╡ 37309807-cd50-4196-a283-473ee937346a
md"## Fit a few absorption lines"

# ╔═╡ 26c1319c-8f39-42d8-a3b7-588ac90054d6
begin
	#λ_lines = [4575.107,4575.542,4575.788, 4576.337,4577.178, 4577.484, 4577.696,4578.032,4578.325]
	λ_lines = [ 4570.0, 4570.8, 4572.4, 4572.6, 4572.8, 4573.0, 4573.2, 4573.5, 4574.1, 4575.5, 4576.0, 4577.1, 4576.35, 4576.8, 4577.6, 4578.4, 4579.8 ]
	σ_lines = fill(0.04,length(λ_lines))
end;

# ╔═╡ 84dcb40f-082f-4867-81fa-67345f96a9c8
md"### Reducing memory allocations"

# ╔═╡ 8fc3ac1f-fe81-46f6-8977-c2291ee8433c
md"### Testing v1"

# ╔═╡ 80403d49-57f3-4b4c-9569-c35244446c1b
md"### Benchmarking implementations"

# ╔═╡ 1d881843-321b-4a96-acac-32ff1ac8bcac
md"""
Obviously, there's room for improvement, but hopefully this demonstrates the basic idea.
You can see that we got pretty good fits when we had a good initial guess for the line location.  For the last line, the initial guess wasn't good and so we have a funky shoulder.  Hopefully, just using line locations from the G2.espresso.mas line list would solve that.  If not, then you could itteratively fit find the bad lines, improve the initial guess and refit.
"""

# ╔═╡ ec7f2e98-b0c8-40ed-a8b1-7d50bbd84503
md"""
# Helper Code
"""

# ╔═╡ fe94411c-4e48-47dc-a358-50a33190cd98
md"## Model for synthetic spectrum"

# ╔═╡ fa1ccd0e-f2b9-4fae-8e2a-205ecf3ce0b4
md"## Blaze function model"

# ╔═╡ 926724fe-e1d9-4409-8b99-ec857e1eb0e5
begin
	struct BlazeModel{T}
		x_mean::T
		polyn::Polynomial{T,:x}
	end
	
	function (blaze::BlazeModel{T})(x::T1)   where { T<:Number, T1<:Number } 
		blaze.polyn(x-blaze.x_mean)
	end
end

# ╔═╡ d537f5a4-50bb-47cc-a960-e9f9b6699ecb
function fit_blaze_model_v0(x::V1, flux::V2, var::V3,
						T::Type = promote_type(T1,T2,T3); order::Integer = 8, mask = 1:length(x)  ) where
			{ T1<:Number, T2<:Number, T3<:Number,
			V1<:AbstractVector{T1}, V2<:AbstractVector{T2}, V3<:AbstractVector{T3} } 
	x_mean = mean(x[mask])
	fit_result = fit( x[mask].-x_mean, 
					  flux[mask], order, weights=1.0./var[mask] )
	return BlazeModel{T}(x_mean,fit_result)
end

# ╔═╡ 8ab5e521-395b-4539-9617-67a8b64011af
function fit_blaze_model(x::V1, flux::V2, var::V3,
						T::Type = promote_type(T1,T2,T3); order::Integer = 8, mask = 1:length(x)  ) where
			{ T1<:Number, T2<:Number, T3<:Number,
			V1<:AbstractVector{T1}, V2<:AbstractVector{T2}, V3<:AbstractVector{T3} } 
	x_mean = mean(view(x,mask))
	fit_result = fit( view(x,mask).-x_mean,
					  view(flux,mask), order, weights=1.0./view(var,mask) )
	return BlazeModel{T}(x_mean,fit_result)
end

# ╔═╡ d25d1008-fd70-4fdb-8adf-8710cc938eb0
# Simplistic fit of polynomial to take out broad pattern to blaze function of spectrograph.  For demo, I just did two itteration of clipping low points (i.e., very likely to be in lines)
# There are much fancier strategies E.g., https://github.com/RvSpectML/NeidSolarScripts.jl/blob/main/src/continuum_rassine_like.jl  But that's likely more detailed than you need to be for this project and trying ot use it could create a distraction.
blaze_model0 =  fit_blaze_model(1:length(λ),flux,var,order=8, mask=mask_fit)

# ╔═╡ 984b85cd-e109-4a00-8c6a-c0efa7dcf35e
begin  # Refit to pixels chosen to not be obviously in lines
	pix_gt_100p_it1 = flux[pix_fit]./blaze_model0.(pix_fit) .>= 1.0
	blaze_model1 =  fit_blaze_model( (1:length(λ)), flux, var,  order=8, mask=pix_fit[pix_gt_100p_it1])
end

# ╔═╡ 7cb4a165-89a2-4020-b74d-ab01ed308965
begin  # Refit to pixels chosen to not be obviously in lines
	pix_gt_100p_it2 = flux[pix_fit]./blaze_model1.(pix_fit) .>= 1.0
	blaze_model2 =  fit_blaze_model( (1:length(λ)), flux, var,  order=8, mask=pix_fit[pix_gt_100p_it2])
end

# ╔═╡ 1baecae4-72cd-43fb-9b27-c619133eb915
begin
	local plt = plot()
	plot!(plt,λ[pix_fit],flux[pix_fit], label="Observation",color=:grey)
	#plot!(λ[pix_fit],blaze_model0.(pix_fit), label="Blaze 0",color=:red)
	#plot!(λ[pix_fit],blaze_model1.(pix_fit), label="Blaze 1",color=:green)
	plot!(λ[pix_fit],blaze_model2.(pix_fit), label="Blaze 2",color=:blue)
	xlabel!("λ (Å)")
	ylabel!("Flux")
end

# ╔═╡ c28ee2aa-d3ee-4688-b9a2-cb6b04f31de8
begin
	local plt = plot(legend=:bottomright)
	#plot!(plt,λ[pix_fit],flux[pix_fit]./blaze_model0.(pix_fit), label="Observation/Blaze 0",color=:red)
	#plot!(plt,λ[pix_fit],flux[pix_fit]./blaze_model1.(pix_fit), label="Observation/Blaze 1",color=:green)
	plot!(plt,λ[pix_fit],flux[pix_fit]./blaze_model2.(pix_fit), label="Observation/Blaze 2", color=:blue)
	xlabel!("λ (Å)")
	ylabel!("Normalized Flux")
end

# ╔═╡ 4e568854-a031-43b8-a2fc-ea5f6da0b89f
# Build dataframe of pixels for ploting (and eventually fitting spectral lines to) in window selected by pix_plt
df = DataFrame( λ=view(λ,pix_plt), 
				flux=view(flux,pix_plt)./blaze_model2.(pix_plt),
				var =view(var,pix_plt)./(blaze_model2.(pix_plt)).^2
				)

# ╔═╡ 92fbbef4-7df3-4f6c-adda-3f8737429f61
begin  # Make workspaces
	dm_tmp = zeros(length(df.λ),1+length(λ_lines)*4)
	covar_tmp = PDiagMat(df.var)
	denom_tmp = zeros(1+length(λ_lines)*4,1+length(λ_lines)*4)
	num_tmp = zeros(1+length(λ_lines)*4)
end

# ╔═╡ 9adb91d3-cffe-4226-b1db-ef100fcbee40
with_terminal() do
	@time fit_blaze_model(1:length(λ),flux,var,order=8, mask=mask_fit)
	@time fit_blaze_model(1:length(λ),flux,var,order=8, mask=mask_fit)
end

# ╔═╡ 890f0d73-3fa4-4ea6-a00a-c3672e8d0fa2
md"""
## Gauss-Hermite parameterization of absorption line
"""

# ╔═╡ 9b3c542e-579b-4f66-8518-4e234cd7c0e7
begin
	struct AbsorptionLine
		λ::Float64
		σ::Float64
		gh_coeff::SVector{4,Float64}  # Statically allocated to reduce memory allocations
	end

	function gaussian(line::AbsorptionLine, λ::Number)
		exp(-((λ-line.λ)/line.σ)^2//2)
	end
	
	# Precompute Gauss-Hermite polynomials once.
	# Should move to module, so not in global scope
	gh_polynomials = [basis(Hermite, i)(variable(Polynomial{Float64})) for i in 0:10]
	
	function gauss_hermite_basis(line::AbsorptionLine, λ, order::Integer)
		@assert 1 <= order+1 <= length(gh_polynomials)
		T = typeof(line.λ)
		gh_poly::Polynomial{T,:x} = gh_polynomials[order+1] 
		x = (λ.-line.λ)./line.σ
		exp.(-0.5.*x.^2) .* gh_poly.(x)		
	end
	
	function gauss_hermite_basis(line::AbsorptionLine, λ; orders::AbstractVector{Int64} = 1:length(line.gh_coeff) )
		T = typeof(line.λ)
		gh_polys::Vector{Polynomial{T,:x}} = gh_polynomials
		x = (λ.-line.λ)./line.σ
		g = exp.(-0.5.*x.^2) 
		p = mapreduce(i->gh_polys[i].(x),hcat,orders)
		g.*p
	end
	
	function (line::AbsorptionLine)(λ)
		T = typeof(line.λ)
		gh_polys::Vector{Polynomial{T,:x}} = gh_polynomials 
		x = (λ.-line.λ)./line.σ
		g = exp.(-0.5.*x.^2) 
		# Could eliminate some allocations in mapreduce with hcat.
		p = mapreduce(i->line.gh_coeff[i].*gh_polys[i].(x),hcat,1:length(line.gh_coeff))
		one(λ) .+ sum(g.*p)
	end

end

# ╔═╡ 8774a77f-bdb6-4ea4-a40c-e96695f7d3e3
function fit_line_v0(λ_line::Number, σ_line::Number, λ::V1, flux::V2, var::V3, T::Type = promote_type(T1,T2,T3); order::Integer=4 ) where
			{ T1<:Number, T2<:Number, T3<:Number,
			  V1<:AbstractVector{T1}, V2<:AbstractVector{T2}, V3<:AbstractVector{T3} } 
	@assert size(λ) == size(flux) == size(var)
	n = length(λ)
	@assert n > 4  # Minimum for 4 parameter fit
	covar = PDiagMat(var)   
	line = AbsorptionLine(λ_line, σ_line, zeros(order) )
	design_matrix = hcat(ones(n), gauss_hermite_basis(line,λ,orders=1:order)) 
	# Generalized linear least squares
	Xt_inv_covar_X = design_matrix' * (covar \ design_matrix) 
	X_inv_covar_y =   design_matrix' * (covar \ flux) 
	coeff_hat =  Xt_inv_covar_X \ X_inv_covar_y
	line = AbsorptionLine(λ_line, σ_line, coeff_hat[2:end] )
	
	#finding the sum of residuals from λ-3σ to λ+3σ
	#need to find the index i_1 in λ where λ[i_1] = λ_line, index i_0 in λ where λ[i_0] = λ_line-3σ, and index i_2 in λ where λ[i_2] = λ_line+3σ
	
	
	return line
end

# ╔═╡ a17f0654-0038-4320-85ca-65221ada0e23
# Try fitting one line at a time
#line = fit_line_v0(4575.5,0.04,df.λ,df.flux,df.var)
#line = fit_line_v0(4576.0,0.04,df.λ,df.flux,df.var)
#line = fit_line_v0(4577.6,0.04,df.λ,df.flux,df.var)
#line = fit_line_v0(4578.46,0.04,df.λ,df.flux,df.var)
line = fit_line_v0(4579.8,0.04,df.λ,df.flux,df.var)


# ╔═╡ cd508e5b-e62d-4c4f-9550-8d9ed3ef2d60
begin
	local plt = plot(legend=:bottomright)
	plot!(plt,df.λ,df.flux, label="Observation/Blaze")
	plot!(plt,df.λ,line.(df.λ),label="Model")
	xlabel!("λ (Å)")
	ylabel!("Normalized Flux")
end

# ╔═╡ a40ffdf7-b2e6-42ee-b4ec-e430fdaed23a
begin
	with_terminal() do
		println(line)
	end
end

# ╔═╡ f5d6e755-140a-44fa-9b9f-6e73227ee9cb
begin
	struct SpectrumModel
		norm::Float64
		lines::Vector{AbsorptionLine}
	end
	
	function (model::SpectrumModel)(λ)
		result = fill(model.norm,length(λ))
		for i in 1:length(model.lines)
			result .*= model.lines[i].(λ)
		end
		return result
	end
end

# ╔═╡ d09a3103-c466-4eb7-8745-7cd1366d6beb
function fit_lines_v0(λ_lines::V1, σ_lines::V2, λ::V3, flux::V4, var::V5, T::Type = promote_type(T1,T2,T3,T4,T5); order::Integer=4 ) where
			{ T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:Number,
			  V1<:AbstractVector{T1}, V2<:AbstractVector{T2}, V3<:AbstractVector{T3}, V4<:AbstractVector{T4}, V5<:AbstractVector{T5}  } 
	@assert size(λ_lines) == size(σ_lines)
	@assert size(λ) == size(flux) == size(var)
	n_pix = length(λ)
	n_lines = length(λ_lines)
	@assert n_lines >= 1 
	@assert 1 <= order <= length(gh_polynomials) 
	@assert n_pix > 1 + order*n_lines  # number of fit parameters
	covar = PDiagMat(var)   # diagonal covariance matrix
	design_matrix = ones(n_pix,1)
	for i in 1:n_lines
		line = AbsorptionLine(λ_lines[i], σ_lines[i],(@SVector zeros(order)) ) #create a line data-structure 
		design_matrix = hcat(  design_matrix,
			gauss_hermite_basis(line,λ,orders=1:order)  )  	# fit to the line	
	end
	Xt_inv_covar_X = design_matrix' * (covar \ design_matrix) 
	X_inv_covar_y =   design_matrix' * (covar \ flux ) 
	coeff_hat = (Xt_inv_covar_X \ X_inv_covar_y)
	
	lines = map(i->AbsorptionLine(λ_lines[i],σ_lines[i], coeff_hat[(2+(i-1)*order):(1+i*order)]), 1:n_lines)
	return SpectrumModel(coeff_hat[1],lines)
end

# ╔═╡ fc6657f6-adf0-4786-8b00-a436e4aa0afe
result0 = fit_lines_v0(λ_lines,σ_lines,df.λ,df.flux,df.var)

# ╔═╡ d35fc6f2-ed06-4b30-bcce-f40d2199c650
with_terminal() do 
	fit_lines_v0, λ_lines,σ_lines,df 
	@time fit_lines_v0(λ_lines,σ_lines,df.λ,df.flux,df.var)
	@time fit_lines_v0(λ_lines,σ_lines,df.λ,df.flux,df.var)
end

# ╔═╡ 17626bf6-08db-413c-9930-28b43b4c1ad9
# Version for writing basis vectors in-place to reduce memory allocations
function gauss_hermite_basis!(output::AA, line::AbsorptionLine, λ::AV2; orders::AVI = 1:length(line.gh_coeff) ) where { T1<:Number, T2<:Number, I<:Integer, AA<:AbstractArray{T1,2}, AV2<:AbstractVector{T2}, AVI<:AbstractVector{I} }
		@assert size(output,1) == length(λ)
		@assert size(output,2) == length(orders)
		T = typeof(line.λ)
		gh_polys::Vector{Polynomial{T,:x}} = gh_polynomials
		x = (λ.-line.λ)./line.σ
		g = exp.(-0.5.*x.^2) 
		p = mapreduce(i->gh_polys[i].(x),hcat,orders)
		output .= g.*p
end
	

# ╔═╡ 763b611a-37a1-41cc-a823-cf384501c03c
function fit_lines_v1(λ_lines::V1, σ_lines::V2, λ::V3, flux::V4, var::V5, T::Type = promote_type(T1,T2,T3,T4,T5); order::Integer=4, covar::AbstractPDMat = PDiagMat(var), prealloc_design_matrix::AA = zeros(T,length(λ),1+length(λ_lines)*order), prealloc_denom::AA = zeros(T,1+length(λ_lines)*order,1+length(λ_lines)*order), prealloc_num::V6 = zeros(T,1+length(λ_lines)*order)) where
			{ T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:Number, T6<:Number,
T7<:Number,	V1<:AbstractVector{T1}, V2<:AbstractVector{T2}, V3<:AbstractVector{T3}, V4<:AbstractVector{T4}, V5<:AbstractVector{T5}, V6<:AbstractVector{T6}, AA<:AbstractArray{T7,2}  } 
	@assert size(λ_lines) == size(σ_lines)
	@assert size(λ) == size(flux) == size(var)
	n_pix = length(λ)
	n_lines = length(λ_lines)
	@assert n_lines >= 1 
	@assert 1 <= order <= length(gh_polynomials) 
	@assert n_pix > 1 + order*n_lines  # number of fit parameters

	design_matrix = prealloc_design_matrix 
	design_matrix[:,1] .= 1
	for i in 1:n_lines
		idx = (2+(i-1)*order):(1+i*order)
		line = AbsorptionLine(λ_lines[i], σ_lines[i],(@SVector zeros(order)) )
		gauss_hermite_basis!(view(design_matrix,:,idx),line,λ,orders=1:order) 
	end
	# Generlized Linear Least Squares
	Xt_inv_covar_X = prealloc_denom
	Xt_inv_covar_X .= design_matrix' * (covar \ design_matrix) 
	X_inv_covar_y = prealloc_num 
	X_inv_covar_y .=   design_matrix' * (covar \ (flux)) 
	coeff_hat =  Xt_inv_covar_X \ X_inv_covar_y
	
	lines = map(i->AbsorptionLine(λ_lines[i],σ_lines[i], view(coeff_hat,(2+(i-1)*order):(1+i*order))), 1:n_lines)
	return SpectrumModel(coeff_hat[1],lines)
end

# ╔═╡ aac63525-a431-4242-b608-b977d2957a55
result = fit_lines_v1(λ_lines,σ_lines,df.λ,df.flux,df.var)

# ╔═╡ d45cc43c-2571-47d0-8ba9-ebbca0235296
begin
	result, result0
	# Comapre to v0
	@test result.norm == result0.norm && 
			all(map(i->result.lines[i] == result0.lines[i],1:length(result.lines)))
end

# ╔═╡ 27b02672-1c6c-449d-ab7f-5b9dd18bd8a1
begin
	local plt = plot(legend=:bottomright)
	plot!(plt,df.λ,df.flux, label="Observation/Fit")
	plot!(plt,df.λ,result(df.λ),label="Model")
	xlabel!("λ (Å)")
	ylabel!("Normalized Flux")
end

# ╔═╡ 6d0abceb-1227-4c58-b50a-80969e451458
result1 = fit_lines_v1(λ_lines,σ_lines,df.λ,df.flux,df.var, covar=covar_tmp, prealloc_design_matrix=dm_tmp, prealloc_denom=denom_tmp,prealloc_num=num_tmp)

# ╔═╡ 36428c0b-bd9e-4f84-9ddf-e52f1324bb5c
begin
	result, result1
	# Compare v1 to itself without preallocating
	@test result.norm == result1.norm &&
      all(map(i->result.lines[i] == result1.lines[i],1:length(result.lines)))
end

# ╔═╡ 38b7abe2-ffed-4af8-a21a-fe4ead027170
with_terminal() do 
	fit_lines_v1,λ_lines,σ_lines,df,covar_tmp, dm_tmp, denom_tmp,num_tmp
	@time fit_lines_v1(λ_lines,σ_lines,df.λ,df.flux,df.var,covar=covar_tmp, prealloc_design_matrix=dm_tmp, prealloc_denom=denom_tmp,prealloc_num=num_tmp)
	@time fit_lines_v1(λ_lines,σ_lines,df.λ,df.flux,df.var,covar=covar_tmp, prealloc_design_matrix=dm_tmp, prealloc_denom=denom_tmp,prealloc_num=num_tmp)
end

# ╔═╡ a0cdca55-c6db-4f9e-90a3-89704637f62b
l_test = AbsorptionLine(4577.6,0.04,(@SVector [1.0,0.5,0.3,0.2]) )

# ╔═╡ 061f8147-a6f0-4a1b-9c19-bd65bc37bcee
md"## Packages used"

# ╔═╡ 4bb36c6c-7be8-4192-8829-e4a12e7adbcd
md"## Other"

# ╔═╡ 99e77b5f-fd69-40cc-b24d-b20d8b348fd9
md"## Code I started on, but ran out of time to make work"

# ╔═╡ f8c5c582-c7ff-4345-975b-0606dc5a6259
function fit_lines_v2(λ_lines::V1, σ_lines::V2, λ::V3, flux::V4, var::V5, T::Type = promote_type(T1,T2,T3,T4,T5); order::Integer=4, covar::AbstractPDMat = PDiagMat(ones(length(λ))), prealloc_design_matrix::AA = zeros(T,length(λ),1+length(λ_lines)*order), prealloc_denom::AA = zeros(T,1+length(λ_lines)*order,1+length(λ_lines)*order), prealloc_num::V6 = zeros(T,1+length(λ_lines)*order)) where
			{ T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:Number, T6<:Number,
T7<:Number,	V1<:AbstractVector{T1}, V2<:AbstractVector{T2}, V3<:AbstractVector{T3}, V4<:AbstractVector{T4}, V5<:AbstractVector{T5}, V6<:AbstractVector{T6}, AA<:AbstractArray{T7,2}  } 
	@assert size(λ_lines) == size(σ_lines)
	@assert size(λ) == size(flux) == size(var)
	n_pix = length(λ)
	n_lines = length(λ_lines)
	@assert n_lines >= 1 
	@assert 1 <= order <= length(gh_polynomials) 
	@assert n_pix > 1 + order*n_lines  # number of fit parameters

	design_matrix = ApplyArray(hcat, ones(n_pix))
	#, randn(3,10))
	#design_matrix = prealloc_design_matrix 
	#design_matrix[:,1] .= 1
	for i in 1:n_lines
		idx = (2+(i-1)*order):(1+i*order)
		line = AbsorptionLine(λ_lines[i], σ_lines[i],(@SVector zeros(order)) )
		design_matrix = ApplyArray(hcat, design_matrix, 
			gauss_hermite_basis(line,λ,orders=1:order) ) 
	end
	#return design_matrix
	# Generlized Linear Least Squares
	Xt_inv_covar_X = prealloc_denom
	Xt_inv_covar_X .= @~ design_matrix' * (covar \ design_matrix) 
	X_inv_covar_y = prealloc_num 
	tmp = covar \ (flux)
	X_inv_covar_y .=    design_matrix' * tmp 
	
	coeff_hat =  collect(Xt_inv_covar_X \ X_inv_covar_y)
	lines = map(i->AbsorptionLine(λ_lines[i],σ_lines[i], view(coeff_hat,(2+(i-1)*order):(1+i*order))), 1:n_lines)
	return SpectrumModel(coeff_hat[1],lines)

end

# ╔═╡ 8e98f9e9-1c28-4566-aeb0-4e1ec38c9467
#result2 = fit_lines_v2(λ_lines,σ_lines,df.λ,df.flux,df.var, covar=covar_tmp, prealloc_design_matrix=dm_tmp, prealloc_denom=denom_tmp,prealloc_num=num_tmp)

# ╔═╡ 114363a2-b696-4046-9b5b-de886520dd27
#@test result2 == result0

# ╔═╡ 7713fa3d-97b9-4be0-894e-1e1eb1c32a00
#=
with_terminal() do 
	λ_lines,σ_lines,df,covar_tmp, dm_tmp, denom_tmp,num_tmp
	@time fit_lines_v2(λ_lines,σ_lines,df.λ,df.flux,df.var,covar=covar_tmp, prealloc_design_matrix=dm_tmp, prealloc_denom=denom_tmp,prealloc_num=num_tmp)
	@time fit_lines_v2(λ_lines,σ_lines,df.λ,df.flux,df.var,covar=covar_tmp, prealloc_design_matrix=dm_tmp, prealloc_denom=denom_tmp,prealloc_num=num_tmp)
end
=#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
FITSIO = "525bcba6-941b-5504-bd06-fd0dc1a4d2eb"
FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
SpecialPolynomials = "a25cea48-d430-424a-8ee7-0d3ad3742e9e"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[compat]
DataFrames = "~1.2.2"
FITSIO = "~0.16.9"
FillArrays = "~0.12.3"
LazyArrays = "~0.21.18"
PDMats = "~0.11.1"
Plots = "~1.21.3"
PlutoTest = "~0.1.0"
PlutoUI = "~0.7.9"
Polynomials = "~2.0.14"
SpecialPolynomials = "~0.2.4"
StaticArrays = "~1.2.12"
StructArrays = "~0.6.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "0f7998147ff3d112fad027c894b6b6bebf867154"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.7.3"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CFITSIO]]
deps = ["CFITSIO_jll"]
git-tree-sha1 = "4379a2dac795014534b9895a45889aa658fca213"
uuid = "3b1b4be9-1499-4b22-8d78-7db3344d1961"
version = "1.4.0"

[[CFITSIO_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Pkg"]
git-tree-sha1 = "2fabb5fc48d185d104ca7ed7444b475705993447"
uuid = "b3e40c51-02ae-5482-8a39-3ace5868dcf4"
version = "3.49.1+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "bdc0937269321858ab2a4f288486cb258b9a0af7"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.3.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "9995eb3977fbf67b86d0a0a0508e83017ded03f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.14.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "727e463cfebd0c7b999bbf3e9e7e16f254b94193"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.34.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "bec2532f8adb82005476c141ec23e921fc20971b"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.8.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d785f42445b63fc86caa08bb9a9351008be9b765"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.2.2"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "fe385ec95ac5533650fb9b1ba7869e9bc28cdd0a"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.5"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FITSIO]]
deps = ["CFITSIO", "Printf", "Reexport", "Tables"]
git-tree-sha1 = "ba5eb4020e474b1c1d4952f91dd7bbfb97b5bf98"
uuid = "525bcba6-941b-5504-bd06-fd0dc1a4d2eb"
version = "0.16.9"

[[FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "5829b25887e53fb6730a9df2ff89ed24baa6abf6"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "0.4.7"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "a3b7b041753094f3b17ffa9d2e2e07d8cace09cd"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.3"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "182da592436e287758ded5be6e32c406de3a2e47"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.58.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "ef49a187604f865f4708c90e3f431890724e9012"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.59.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "60ed5f1643927479f845b0135bb369b031b541fa"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.14"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "8a954fed8ac097d5be04921d595f741115c1b2ad"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+0"

[[HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "81e1680c7242061bca2c8e5614d36802a97d1999"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.5"

[[HypertextLiteral]]
git-tree-sha1 = "1e3ccdc7a6f7b577623028e0095479f4727d8ec1"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.8.0"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IrrationalConstants]]
git-tree-sha1 = "f76424439413893a832026ca355fe273e93bce94"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "MatrixFactorizations", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "4a67e454dbcf6677e2ab7a0de318ea9e6f21169a"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "0.21.18"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "3d682c07e6dd250ed082f883dc88aee7996bf2cc"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.0"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "0fb723cd8c45858c22169b2e42269e53271a6df7"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.7"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "24814f4e65b4521ba081ccaaea9f5c6533c462a2"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "0.8.4"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Memoize]]
deps = ["MacroTools"]
git-tree-sha1 = "2b1dfcba103de714d31c033b5dacc2e4a12c7caa"
uuid = "c03570c3-d221-55d1-a50c-7939bbd78826"
version = "0.4.4"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "2ca267b08821e86c5ef4376cffed98a46c2cb205"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "748f6e1e4de814b101911e64cc12d83a6af66782"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.2"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9ff1c70190c1c30aebca35dc489f7411b256cd23"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.13"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "2dbafeadadcf7dadff20cd60046bba416b4912be"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.21.3"

[[PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "3479836b31a31c29a7bac1f09d95f9c843ce1ade"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.1.0"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "0bbfdcd8cda81b8144de4be8a67f5717e959a005"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.14"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a193d6ad9c45ada72c14b731a318bedd3c2f00cf"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.3.0"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "0d1245a357cc61c8cd61934c07447aa569ff22e6"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.1.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "d4491becdc53580c6dadb0f6249f90caae888554"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.0"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "a322a9493e49c5f3a10b50df3aedaf1cdb3244b7"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.1"

[[SpecialPolynomials]]
deps = ["FastGaussQuadrature", "HypergeometricFunctions", "Intervals", "LinearAlgebra", "Memoize", "Polynomials", "QuadGK", "Requires", "SpecialFunctions"]
git-tree-sha1 = "de22f4e1699eb069eb448b5970675ae2049e41cd"
uuid = "a25cea48-d430-424a-8ee7-0d3ad3742e9e"
version = "0.2.4"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8cbbc098554648c84f79a463c9ff0fd277144b6c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.10"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "1700b86ad59348c0f9f68ddc95117071f947072d"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.1"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "d0c690d37c73aeb5ca063056283fde5585a41710"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "6c9040665b2da00d30143261aea22c7427aada1c"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.7"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─39311305-5492-441b-8346-23d783a04220
# ╠═e624fdcd-df68-4cb7-811a-80154f624a47
# ╠═93b5f3b9-8d7a-43e4-9901-d0f8800b2b28
# ╠═79f795c4-a1a6-4e08-a3f1-3bd7bca450ba
# ╠═0a132ecb-c759-4448-823a-9c4727ee6916
# ╠═1870ae31-79bb-4e5d-a7f1-6202cca94867
# ╟─1d610ca0-6981-48ff-b9d8-76c45807b5e7
# ╠═d25d1008-fd70-4fdb-8adf-8710cc938eb0
# ╠═1baecae4-72cd-43fb-9b27-c619133eb915
# ╠═984b85cd-e109-4a00-8c6a-c0efa7dcf35e
# ╠═7cb4a165-89a2-4020-b74d-ab01ed308965
# ╠═c28ee2aa-d3ee-4688-b9a2-cb6b04f31de8
# ╠═c9a3550e-ab84-44d4-a935-64e80ed51d63
# ╠═4e568854-a031-43b8-a2fc-ea5f6da0b89f
# ╟─8d935f79-f2da-4d41-8dad-85cd08197d17
# ╠═8774a77f-bdb6-4ea4-a40c-e96695f7d3e3
# ╠═a17f0654-0038-4320-85ca-65221ada0e23
# ╠═cd508e5b-e62d-4c4f-9550-8d9ed3ef2d60
# ╠═a40ffdf7-b2e6-42ee-b4ec-e430fdaed23a
# ╟─37309807-cd50-4196-a283-473ee937346a
# ╠═26c1319c-8f39-42d8-a3b7-588ac90054d6
# ╠═d09a3103-c466-4eb7-8745-7cd1366d6beb
# ╟─84dcb40f-082f-4867-81fa-67345f96a9c8
# ╠═763b611a-37a1-41cc-a823-cf384501c03c
# ╟─8fc3ac1f-fe81-46f6-8977-c2291ee8433c
# ╠═aac63525-a431-4242-b608-b977d2957a55
# ╠═fc6657f6-adf0-4786-8b00-a436e4aa0afe
# ╠═d45cc43c-2571-47d0-8ba9-ebbca0235296
# ╠═92fbbef4-7df3-4f6c-adda-3f8737429f61
# ╠═6d0abceb-1227-4c58-b50a-80969e451458
# ╠═36428c0b-bd9e-4f84-9ddf-e52f1324bb5c
# ╟─80403d49-57f3-4b4c-9569-c35244446c1b
# ╠═d35fc6f2-ed06-4b30-bcce-f40d2199c650
# ╠═38b7abe2-ffed-4af8-a21a-fe4ead027170
# ╠═27b02672-1c6c-449d-ab7f-5b9dd18bd8a1
# ╟─1d881843-321b-4a96-acac-32ff1ac8bcac
# ╠═ec7f2e98-b0c8-40ed-a8b1-7d50bbd84503
# ╠═fe94411c-4e48-47dc-a358-50a33190cd98
# ╠═f5d6e755-140a-44fa-9b9f-6e73227ee9cb
# ╟─fa1ccd0e-f2b9-4fae-8e2a-205ecf3ce0b4
# ╠═926724fe-e1d9-4409-8b99-ec857e1eb0e5
# ╠═d537f5a4-50bb-47cc-a960-e9f9b6699ecb
# ╠═8ab5e521-395b-4539-9617-67a8b64011af
# ╠═9adb91d3-cffe-4226-b1db-ef100fcbee40
# ╟─890f0d73-3fa4-4ea6-a00a-c3672e8d0fa2
# ╠═9b3c542e-579b-4f66-8518-4e234cd7c0e7
# ╠═17626bf6-08db-413c-9930-28b43b4c1ad9
# ╠═a0cdca55-c6db-4f9e-90a3-89704637f62b
# ╟─061f8147-a6f0-4a1b-9c19-bd65bc37bcee
# ╠═caef09cc-0e00-11ec-1753-d7e117eb8c20
# ╠═10393148-b9b0-44a0-9b47-c6780318316b
# ╠═81ab1954-a65e-4633-a340-2e707e5ce879
# ╠═6b4745bf-e13f-4f20-9bda-4dd37662325b
# ╠═4d1ad1b6-788d-4607-b4c3-b37a4cd37e3e
# ╠═05df8c3f-f679-489f-9e51-e5fad58d9a55
# ╠═83cd8a9e-7bdc-4834-b910-8068767ebfac
# ╠═e533ddde-2ed7-42b3-a1cd-a1e5738cf3b3
# ╠═ce588b8c-4656-4b58-b4e1-9c0ae8b9eefd
# ╟─4bb36c6c-7be8-4192-8829-e4a12e7adbcd
# ╟─99e77b5f-fd69-40cc-b24d-b20d8b348fd9
# ╠═f8c5c582-c7ff-4345-975b-0606dc5a6259
# ╠═8e98f9e9-1c28-4566-aeb0-4e1ec38c9467
# ╠═114363a2-b696-4046-9b5b-de886520dd27
# ╠═7713fa3d-97b9-4be0-894e-1e1eb1c32a00
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
