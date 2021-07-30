# Slice sampler with 0 lb, step size 1, 100 step limit
# Based on http://www.cs.toronto.edu/~radford/slice.software.html
function sliceSample(x::Float32,f)
	y = f(x)-rand(Exponential(1.0))

	left = x-rand()
	right = left+1.0f0
	t = rand()*100
	leftc = floor(t)
	rightc = 99-leftc

	while true
		if left<0.000001f0
			left = 0.000001f0
			break
		end
		if f(left)<=y
			break
		end
		left = left-1.0f0
		leftc = leftc-1
		if leftc<=0
			break
		end
	end
	if left<0.000001f0
		left = 0.000001f0
	end

	while true
		if f(right)<=y
			break
		end
		right = right+1.0f0
		rightc = rightc-1
		if rightc<=0
			break
		end
	end

	xx::Float32 = left+rand()*(right-left)
	fxx = f(xx)
	step = 0
	while step<100 
		if isfinite(fxx)
			if fxx>=y
				break
			end
			if xx>x
				right = xx
			else
				left = xx
			end
		end
		xx = left+rand()*(right-left)
		fxx = f(xx)
		step = step+1
	end
	if !isfinite(fxx) # valid f(x)?
		xx = x
	end
	return xx
end

