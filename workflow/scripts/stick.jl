mutable struct MixtureTheta{T}
	alpha :: Float32
	z :: Array{Int64}
	K :: Int64
	cluster :: Array{T}
	counts :: Array{Int64}
end

struct MixtureDistribution
	f
	priorSample
	priorSample!
	priorDens
	condLocalSample!
	condGlobalSample!
	L :: Int64
	est :: Int64
end

function runif(a,b)
	x = Float32(rand(Uniform(a,b)))
	if x==a
		x = x+floatmin(Float32)
	else
		if x==b
			x = x-floatmin(Float32)
		end
	end
	return x
end

function rbeta(a,b)
	x = Float32(rand(Beta(a,b)))
	if x==0.0
		x = x+floatmin(Float32)
	else
		if x==1.0
			x = x-floatmin(Float32)
		end
	end
	return x
end

function sampleAlpha!(thetaDP::MixtureTheta{T}) where T
	# non-empty	
	k = 0
	for j = 1:thetaDP.K
		ix = findall(thetaDP.z .== j)
		if length(ix)>0
			k = k + 1
		end
	end
	# https://stats.stackexchange.com/questions/381867/gibbs-sampler-for-dirichlet-process-concentration-parameter
	n = size(thetaDP.z)[1]
	ua = n/(n+thetaDP.alpha)
	t = rand()
	if t<ua
		u = 1.0f0
	else
		u = 0.0f0
	end
	v = rand(Beta(thetaDP.alpha+1,n))
	# julia distributions gamma is shape/scale
	alpha = rand(Gamma(1+k-1.0f0+u,1.0f0/(1.0f0-log(v))))
	thetaDP.alpha = alpha
	#println("Alpha "*string(alpha)*", "*string(k)*" clusters")
end

function gibbsSampleDP!(mix,thetaDP::MixtureTheta{T},thetaModel) where T
	
	# STEP A
	# Z* = thetaDP.K
	ZZ = thetaDP.K
	nc = thetaDP.counts

	# STEP B
	# B1
	# Sample V^A
	vA = zeros(Float32,ZZ)
	for c in 1:ZZ
		vA[c] = rbeta(1+nc[c],thetaDP.alpha+sum(nc[c+1:ZZ]))
	end
	# B2
	# Gibbs sampling of active theta - already in outer loop
	# B3
	if ZZ>1
		# Swap move 1
		nonempty = []
		for c in 1:ZZ
			if nc[c]>0
				push!(nonempty,c)
			end
		end
		if length(nonempty)>1
			sw = sample(nonempty,2,replace=false)
			a = sw[1]
			b = sw[2]
			pa = vA[a]*prod(1.0f0.-vA[1:(a-1)])
			pb = vA[b]*prod(1.0f0.-vA[1:(b-1)])
			acc = (nc[b]-nc[a])*(log(pa)-log(pb))
			t = rand()
			if t<exp(acc)
				# swap a and b				
				for i in 1:mix.L
					if thetaDP.z[i]==a
						thetaDP.z[i] = b
					else
						if thetaDP.z[i]==b
							thetaDP.z[i] = a
						end
					end
				end
				tmp__ = thetaDP.cluster[a]
				thetaDP.cluster[a] = thetaDP.cluster[b] # ??
				thetaDP.cluster[b] = tmp__
				tmpn = nc[a]
				nc[a] = nc[b]
				nc[b] = tmpn
			end
		end
		# Swap move 2
		s = sample(1:(ZZ-1))
		acc = nc[s]*log(1.0f0-vA[s+1])-nc[s+1]*log(1.0f0-vA[s])
		t = rand()
		if t<exp(acc)
			# swap s and s+1
			v1 = vA[s]
			v2 = vA[s+1]
			vA[s+1] = v1
			vA[s] = v2
			for i in 1:mix.L
				if thetaDP.z[i]==s
					thetaDP.z[i] = s+1
				else
					if thetaDP.z[i]==(s+1)
						thetaDP.z[i] = s
					end
				end
			end
			tmp__ = thetaDP.cluster[s]
			thetaDP.cluster[s] = thetaDP.cluster[s+1]
			thetaDP.cluster[s+1] = tmp__
			tmpn = nc[s]
			nc[s] = nc[s+1]
			nc[s+1] = tmpn
			if s==(ZZ-1)
				if nc[ZZ]==0
					ZZ = s
				end
			end
		end
	end
	# B4
	# Sample U_i ~ Uniform[0,V_{Z_i}\prod_{l<Z_l} (1-V_l)]
	u = zeros(Float32,mix.L)
	for i in 1:mix.L
		l = thetaDP.z[i]
		vv = vA[l]*prod(1.0f0.-vA[1:(l-1)])
		u[i] = runif(0,vv) #rand(Uniform(0,vv))
	end

	# STEP C
	# Calculate Z*, U*
	uu = minimum(u)

	# STEP D
	# D1
	# Sample DP alpha (skip for now)
	# D2
	# Sample V^P (potential) - here CC is C*, total number of potential + active
	vvsum = vA[1]
	vv = 1.0f0
	for i in 2:ZZ
		vv = vv*(1.0f0-vA[i-1])
		vvsum = vvsum + vA[i]*vv
	end
	CC = ZZ
	vv = vv*(1-vA[ZZ])

	while true 
		if CC>500
			break
		end
		if vvsum>(1-uu)
			break
		end
		if vv==0.0
			break
		end
		CC = CC+1
		push!(vA,rbeta(1,thetaDP.alpha))
		vvsum = vvsum + vA[CC]*vv
		vv = vv*(1.0f0-vA[CC])
	end
	
	# STEP E
	# Sample theta_P
	# delete old...
	thetaDP.cluster = thetaDP.cluster[1:ZZ]
	for j = (ZZ+1):CC
		push!(thetaDP.cluster,mix.priorSample())
	end

	# STEP F
	# Update global parameters - already in outer loop
	
	# STEP G 
	# Sample from vA, vP
	countsnew = zeros(Int64,CC)
	r = zeros(Float32,CC)
	for i = 1:mix.L
		vv = 1.0f0
		for j = 1:CC
			r[j] = log(vA[j]*vv)+mix.f(i,thetaDP.cluster[j],thetaModel)
			vv = vv*(1.0f0-vA[j])
			if isnan(r[j])
				r[j]=-Inf
			end
		end
		
		rr = r.-maximum(r)
		rr = exp.(rr)
		rr = rr./sum(rr)
		thetaDP.z[i] = rand(Categorical(rr))
		countsnew[thetaDP.z[i]] = countsnew[thetaDP.z[i]]+1
	end

	# Update values
	thetaDP.K = maximum(thetaDP.z)
	thetaDP.counts = countsnew
end

function initialiseDP(mix,ninit)
	counts = zeros(Int64,ninit)
	initialZ = zeros(Int64,mix.L)
	for i in 1:mix.L
		initialZ[i] = rand(Categorical(ones(Float32,ninit)/ninit))
		counts[initialZ[i]] = counts[initialZ[i]]+1
	end
	clusters = Array{Any}(undef,ninit)
	for j in 1:ninit
		clusters[j] = mix.priorSample()
	end
	MixtureTheta(1.0f0,initialZ,ninit,clusters,counts)
end

function calcpost(thetaDP,thetaModel,mix)
	p = 0.0f0
	for j = 1:thetaDP.K
		ix = findall(thetaDP.z .== j)
		if length(ix)>0
			for i in ix
				p = p + mix.f(i,thetaDP.cluster[j],thetaModel)
			end
			p = p + mix.priorDens(thetaDP.cluster[j])
			p = p + llgamma(length(ix))
		end
	end
	return p
end

function mcmcDP(mix,thetaModel::M,nsteps,burn,thin,ngibbs,ninit) where {M}
	thetaDP = initialiseDP(mix,ninit)
	sampc = 1
	nsamp = trunc(Int,(nsteps-burn)/thin)
	sampsDP = Array{typeof(thetaDP)}(undef,nsamp)
	sampsModel = Array{M}(undef,nsamp)
	post = Array{Float32}(undef,nsamp)
	for step = 1:nsteps
		# Update cluster parameters
		for j = 1:thetaDP.K
			ix = findall(thetaDP.z .== j)
			if length(ix)>0
				for gsamp in 1:ngibbs
					mix.condLocalSample!(thetaDP.cluster[j],thetaModel,ix)
				end
			else
				mix.priorSample!(thetaDP.cluster[j])
			end
		end
		# sample global model parameters
		mix.condGlobalSample!(thetaDP,thetaModel)
		# sample DP mixture
		gibbsSampleDP!(mix,thetaDP,thetaModel)
		sampleAlpha!(thetaDP)
		if step>burn && ((step-burn)%thin)==0
			sampsDP[sampc] = deepcopy(thetaDP)
			sampsModel[sampc] = deepcopy(thetaModel)
			post[sampc] = calcpost(thetaDP,thetaModel,mix)
			sampc = sampc + 1
		end
	end
	return (sampsDP,sampsModel,post)
end
