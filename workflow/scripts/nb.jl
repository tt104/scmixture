function llgamma(x::Int64)
	return (logabsgamma(x))[1]
end

function llgamma(x::Float32)
	return (logabsgamma(x))[1]
end

function llgamma(x::Float64)
	return (logabsgamma(x))[1]
end

function logsum(x,y)
	m = max(x,y)
	s = log(exp(x-m)+exp(y-m))+m
	return s
end

# Data structures

struct Hyperparameters
	a :: Float32
	b :: Float32
	s :: Float32
	t :: Float32
end

mutable struct ThetaClust
	mu :: Array{Float32,1}
	omega :: Array{Float32,1}
	w :: Array{Float32,1}
end

mutable struct Model
	scales :: Array{Float32,1}
	ziflag :: Int32
end

# Prior

function mkPrior(D,h)
	M = size(D)[2]
	function sampleClusterPrior()
		mu = randGamma(h.s,1.0f0/h.t,M) # Hyperparameters are shape/rate
		omega = randGamma(h.a,1.0f0/h.b,M)
		w = rand(Float32,M)
		ThetaClust(mu,omega,w)
	end
end

function mkPriorDens(h)
	function priorDens(cl)
		l = sum(logpdf.(Gamma(h.s,1.0f0/h.t),cl.mu))
		l = l+sum(logpdf.(Gamma(h.a,1.0f0/h.b),cl.omega))
		return l
	end
end

function mkPrior!(D,h)
	M = size(D)[2]
	function sampleClusterPrior(cl)
		randGamma!(h.s,1.0f0/h.t,cl.mu) # Hyperparameters are shape/rate
		randGamma!(h.a,1.0f0/h.b,cl.omega)
		rand!(cl.w)
	end
end

# Likelihood
function mkLikelihood(D::Array{Float32,2})
	M = size(D)[2]
	function f(i::Int64,theta::ThetaClust,model::Model)
		s = 0.0f0
		delta = x -> x==0.0f0 ? 1.0 : 0.0
		@inbounds for j = 1:M
			io = theta.omega[j]
			mo = model.scales[i]*theta.mu[j]*(1.0f0/theta.omega[j])
			w = theta.w[j]
			l = llgamma(D[i,j]+io) - (llgamma(D[i,j]+1.0f0) + llgamma(io)) - (io)*log(1.0f0+mo) + D[i,j] * (log(mo) - log(1.0f0+mo))
			s = s + logsum(log(w)+l,log(1-w)+log(delta(D[i,j])))
		end
		return s
	end
end

# Local Gibbs update- TO FIX --------------------------------------------------------------------------------------------------------------------------- #

function mkCondLocal(D,h)
	M = size(D)[2]
	function sampleCondLocal(theta,thetaGlobal,ix)
		xs = D[ix,:]
		L = length(ix)
		scalesix = thetaGlobal.scales[ix]
		#sample t
		t = zeros(L,M)
		for i in 1:L
			for j in 1:M
				r0 = delta(xs[i,j])
				m = theta.mu[j]
				io = theta.omega[j]
				mo = scalesix[i]*m*(1.0f0/io)
				l = llgamma(xs[i,j]+io) - (llgamma(xs[i,j]+1.0f0) + llgamma(io)) - (io)*log(1.0f0+mo) + xs[i,j] * (log(mo) - log(1.0f0+mo)) 
				rN = exp(l)
				if r0>0.00
					r = (theta.w[j]*rN)/((theta.w[j]*rN)+((1.0f0-theta.w[j])*r0))
				else
					r = 1.0
				end
				tmp = rand()
				if tmp<r
					t[i,j] = 1.0f0
				else
					t[i,j] = 0.0f0
				end
			end
		end
		#sample v
		v = zeros(L,M)
		for i in 1:L
			for j in 1:M
				v[i,j] = randGamma(xs[i,j]*t[i,j]+theta.omega[j],1.0f0/(scalesix[i]*theta.mu[j]*t[i,j]+theta.omega[j])) # Shape/Rate param
			end
		end
		#sample mu
		for j in 1:M
			theta.mu[j] = randGamma(dot(xs[:,j],t[:,j])+h.s,1.0f0/(dot(scalesix.*v[:,j],t[:,j])+h.t)) # Shape/Rate params
		end
		#sample omega
		for j in 1:M
			N = sum(t[:,j])
			lvSum = 0.0f0
			for i in 1:L
				lvSum += t[i,j]*log(v[i,j])
			end
			vSum = dot(t[:,j],v[:,j])
			function likeOmega(w)::Float32
				return (w-1)*lvSum-w*vSum+(h.a-1.0f0)*log(w)-h.b*w+N*w*log(w)-N*llgamma(w)
			end
			theta.omega[j] = sliceSample(theta.omega[j],likeOmega)
		end
		#sample w
		for j in 1:M
			theta.w[j] = rbeta(1+sum(t[:,j]),1+(M-sum(t[:,j]))) # beta
		end
	end
end

# Global Gibbs update
# empty
function mkCondGlobal(D::Array{Float32,2})
	L,M = size(D)
	function sampleCondGlobal(thetaDP,thetaModel::Model)
	end
end
