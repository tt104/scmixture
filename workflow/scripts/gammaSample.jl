using GSL

# Workaround for https://github.com/JuliaStats/Distributions.jl/issues/1024

GSL_T = GSL.rng_env_setup()
GSL_r = GSL.rng_alloc(GSL_T)

function randGamma(a,b)
	global GSL_r
	x = Float32(GSL.ran_gamma(GSL_r,a,b))
	return x
end

function randGamma(a,b,n)
	global GSL_r
	xs = zeros(Float32,n)
	for i in 1:n
		xs[i] = Float32(GSL.ran_gamma(GSL_r,a,b))
	end
	return xs
end

function randGamma!(a,b,xs)
	global GSL_r
	for i in 1:length(xs)
		xs[i] = Float32(GSL.ran_gamma(GSL_r,a,b))
	end
end
