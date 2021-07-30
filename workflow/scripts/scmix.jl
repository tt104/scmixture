import Pkg
Pkg.activate(snakemake.input[4])
Pkg.instantiate()

using Distributions #
using StatsBase
using Random
using DelimitedFiles #
using SpecialFunctions
using Statistics
using Clustering #
using LinearAlgebra #

include(snakemake.scriptdir*"/sliceSample.jl")
include(snakemake.scriptdir*"/stick.jl")
include(snakemake.scriptdir*"/mult.jl")
include(snakemake.scriptdir*"/gammaSample.jl")
include(snakemake.scriptdir*"/nb.jl")

# Read input data

input_file = snakemake.input[1]
scales_file = snakemake.input[2]
hyper_params = snakemake.input[3]

filter = snakemake.config["filter"]
prior_samps = 1
gibbs_samps = 1
model_type = snakemake.config["model"]

(genes_expr,cell_labels) = readdlm(input_file,',',header=true)

gene_names = genes_expr[:,1]

expr = genes_expr[:,2:(size(genes_expr)[2])]

dataFull = copy(expr')
dataFull = Float32.(dataFull)

println("Reading scales from "*scales_file)
scales = readdlm(scales_file,',',Float32)[:,1]

# Filter by CoV

datascaled = copy(dataFull)
for i in 1:size(datascaled)[1]
	datascaled[i,:] = datascaled[i,:]/scales[i]
end

P = size(datascaled)[2]
vs = zeros(Float64,P)
ldata = log.(datascaled.+1)
for i in 1:P
    vs[i] = var(ldata[:,i])
end
cv = sqrt.(exp.(vs.-1))
ixc = sortperm(cv,rev=true)
data = dataFull[:,ixc[1:filter]]

# K means for initial Z
kmeanN = snakemake.config["kmeans"]
datak = log.(datascaled[:,ixc[1:filter]]'.+1)
kmc = kmeans(datak,kmeanN)
initZ = assignments(kmc)

cells = string(size(data)[1])
transcripts_full = string(size(dataFull)[2])
transcripts = string(size(data)[2])

println("Data with "*cells*" cells and "*transcripts_full*" transcripts, filtered to "*transcripts)

# Output helper functions
function getZs(samps,map_ix)
	zs = samps[map_ix].z
	return zs
end

function getAllZs(samps)
	zs = zeros(length(samps),length(samps[1].z))
	for i in 1:length(samps)
		zs[i,:] = samps[i].z
	end
	return zs
end

function getMus(samps,map_ix)
	i = map_ix
	zs = zeros(length(samps[i].cluster),length(samps[i].cluster[1].mu))
	for j in 1:length(samps[i].cluster) 
		zs[j,:] = samps[i].cluster[j].mu
	end
	return zs
end

function getOmegas(samps,mapix)
	i = map_ix
	zs = zeros(length(samps[i].cluster),length(samps[i].cluster[1].omega))
	for j in 1:length(samps[i].cluster) 
		zs[j,:] = samps[i].cluster[j].omega
	end
	return zs
end

L = size(data)[1]


if model_type=="nb"
	println("Using DP NB model")
	println("Reading hyperparameters from "*hyper_params)
	hp = readdlm(hyper_params,',',Float32)
	h = Hyperparameters(hp[1],hp[2],hp[3],hp[4]) #EB TEST
	println("Using hyperparameters "*string(h))
	mixdist = MixtureDistribution(mkLikelihood(data),mkPrior(data,h),mkPrior!(data,h),mkPriorDens(h),mkCondLocal(data,h),mkCondGlobal(data),L,prior_samps)
	glob = Model(ones(size(data)),scales,1)
else
	println("Using DP multinomial model")
	h = HyperparametersMult(1.0f0)
	mixdist = MixtureDistribution(mkLikelihoodMult(data),mkPriorMult(data,h),mkPriorMult!(data,h),mkPriorDensMult(data,h),mkCondLocalMult(data,h),mkCondGlobalMult(data),L,prior_samps)
	glob = ModelMult(ones(size(data)),scales)
end

burnin = snakemake.config["burnin"]
iter = snakemake.config["iter"]
thin = snakemake.config["thin"]
chains = snakemake.config["chains"]

#Run sampler
map_ix = 0
map_est = -Inf
zmat = zeros(Float32,(L,L))
sample_count = 0
for sampruns in 1:chains
	s = mcmcDP(mixdist,glob,iter,burnin,thin,gibbs_samps,initZ)
	m_ix = argmax(s[3])
	m = maximum(s[3])
	if m>map_est
		global map_est = m
		global map_ix = m_ix
		global samp = s
	end
	# Build matrix
	zs = getAllZs(s[1])
	ns = size(zs)[1]
	global sample_count = sample_count + ns
	for i in 1:(L-1)
		for j in (i+1):L
			for m in 1:ns
				if zs[m,i]==zs[m,j]
					global zmat[i,j] = zmat[i,j]+1
				end
			end
		end
	end
end

writedlm(snakemake.output[1],getZs(samp[1],map_ix),',')

writedlm(snakemake.output[4],gene_names[ixc[1:filter]],',')

if model_type=="nb"
	writedlm(snakemake.output[2],getMus(samp[1],map_ix),',')
	writedlm(snakemake.output[3],getOmegas(samp[1],map_ix),',')
end

# Uncertain labels

zmat = zmat+zmat'
zmat = zmat./sample_count
z = getZs(samp[1],map_ix)

function getUncert(z,zmat)
	L = size(z)[1]
	u = []
	for i in 1:L
		c = z[i]
		ne = findall(z.==c)
		filter!(x->x!=i,ne)
		if length(ne)>0
			ps = zmat[i,ne]
			p = sum(ps)/length(ne)
			if p<0.5
				push!(u,i)
			end
		end
	end
	return u
end

u = getUncert(z,zmat)

writedlm(snakemake.output[5],cell_labels[u],',')
