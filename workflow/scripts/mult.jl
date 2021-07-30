# Data structures

struct HyperparametersMult
	alp :: Float32
end

mutable struct ThetaClustMult
	ps :: Array{Float32,1}
end

mutable struct ModelMult
	U :: Array{Float32,2}
	scales :: Array{Float32,1}
end

# Prior

function mkPriorMult(D,h)
	M = size(D)[2]
	function sampleClusterPrior()
		ps = rand(Dirichlet(M,h.alp))
		ThetaClustMult(ps)
	end
end

function mkPriorDensMult(D,h)
	M = size(D)[2]
	function priorDens(cl)
		l = logpdf(Dirichlet(M,h.alp),cl.ps)
		l
	end
end

function mkPriorMult!(D,h)
	M = size(D)[2]
	function sampleClusterPrior(cl)
		rand!(Dirichlet(M,h.alp),cl.ps)
	end
end

# Likelihood
function mkLikelihoodMult(D::Array{Float32,2})
	M = size(D)[2]
	function f(i::Int64,theta::ThetaClustMult,model::ModelMult)
		n = Int(sum(D[i,:]))
		s = logpdf(Multinomial(n,theta.ps),D[i,:])
		return s
	end
end

# Local Gibbs update
function mkCondLocalMult(D,h)
	M = size(D)[2]
	function sampleCondLocal(theta,thetaGlobal,ix)
		xs = D[ix,:]
		L = length(ix)
		rand!(Dirichlet(h.alp.+vec(sum(xs,dims=1))),theta.ps)
	end
end

# Global Gibbs update
function mkCondGlobalMult(D::Array{Float32,2})
	L,M = size(D)
	function sampleCondGlobal(thetaDP,thetaModel::ModelMult)
	end
end
