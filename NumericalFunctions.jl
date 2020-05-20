using Polynomials
using Plots
using LinearAlgebra

##---------- Fixed Point Iteration ----------##

# Takes: fixed point function, the point of interest and # of iterations
function fpi(g::Function, x::Float64, n::Int64) 

	for j = 1:n
		gx = g(x) # evaluate the function at desired point
		x = gx # update the old point
	end
	return ["j" "x";j x] # return the convergence of the point
end

##---------- Bisection Method ----------##

# Takes: function, range [a, b], and the level of tolerance for error
function bisec(f::Function, a::Int8, b::Int8, tol::Float64) 
	
	# Compute Initial function values
	fa = f(a)
	fb = f(b)
	fmid = 0

	if(fa * fb < 0) # function values must be opposite signs for this to work
		for j = 1:10000
			
			# compute midpoint of a & b & calculate function at mid
			mid = a + (b - a) / 2 
			fmid = f(mid)
			
			# if the function at mid is 0 or the distance from a to b is below
			# tolerance level, then you've found the point & end
			if (fmid == 0) || ((b - a) / 2 < tol)
				return mid
			end
			
			# if function at mid has the same sign as function at a,
			# use mid as the new a, otherwise use mid as the new b
			if fmid * fa > 0
				a = mid
				fa = fmid
			else
				b = mid
			end
		end
	end
end

##---------- Linear Equation Solver ----------##

# small simple solver for linear equations
function LES(a::Float64, b::Float64, e::Float64, c::Float64, d::Float64, f:Float64)

	m = c / a
	d1 = d - m * b
	f1 = f - m * e
	y = f1 / d1
	x = (e - b * y) / a
	return [x y]
end

##---------- Combinations ----------##

# recursively creates all the combinations of the elements of x
function combos(x::Array{Float64, 1}, c::Array{Float64, 1}, i::Int8, j::Int8)
	
		if i == length(x)-1
		
			return push!(c,x[i]*x[j])
		
		elseif j < length(x)
	
			perm(x, push!(c,x[i]*x[j]), i, j+=1)
		
		elseif j == length(x)
		
			perm(x, push!(c,x[i]*x[j]), i+=1, i+1)
		end
end

##---------- Find Min Index ----------##

# two functions that find the minimum value in a given function and returns the x value at which it occurs. One function is the absolute valued min.
# Takes: a start and an end or a vector of points, and a function
function mindex(a::Float64, b::Float64, f::Function)
	
	min = f(a)
	dex = a

	for x = LinRange(a, b, 1000)
		if min > f(x)
			min = f(x)
			dex = x
		end
	end
	return [min; dex]
end

function mindex(x::Array{Float64, 1}, f::Function)

	min = abs.(f(x[1]))
	dex = x[1]

	for j = 2:length(x)
		if min > abs.(f(x[j]))
			min = abs.(f(x[j]))
			dex = x[j]
		end
	end
	return [min; dex]
end

##---------- Find Max Index ----------##

# same as the find min index, but for the max....
function maxdex(a::Float64, b::Float64, f::Function)

	max = f(a)
	dex = a

	for x = LinRange(a, b, 1000)
		if max < f(x)
			max = f(x)
			dex = x
		end
	end
	return [max; dex]
end

function maxdex(x::Array{Float64, 1}, f::Function)

	max = abs.(f(x[1]))
	dex = x[1]

	for j = 2:length(x)
		if max < abs.(f(x[j]))
			max = abs.(f(x[j]))
			dex = x[j]
		end
	end
	return [max; dex]
end

##---------- Euler's Method ----------##

function eulr(dy::Function, t::Array{Float64, 1}, y0::Float64, y_eulr::Array{Float64, 1})

	y_eulr[1] = y0
	
	for j = 1:length(t) - 1

		y_eulr[j + 1] = y_eulr[j] + h*dy(t[j], y_eulr[j]) # new = old + roc*timestep
	end
	
	return y_eulr
end

##---------- Midpoint Method ----------##

function midpoint(dy::Function, t::Array{Float64, 1}, y0::Float64, y_mid::Array{Float64, 1})
	
	y_mid[1] = y0
	
	for j = 1:length(t) - 1
		
		roc = dy(t[j] + h/2, y_mid[j] + (h/2)*dy(t[j], y_mid[j])) # half an Euler's step
		y_mid[j + 1] = y_mid[j] + h*roc

	end

	return y_mid
end

##---------- Runge Kutta (4th Order) ----------##

function rk4(dy::Function, t::Array{Float64, 1}, y0::Float64, y_rk4::Array{Float64, 1})

	y_rk4[1] = y0

	for j = 1:length(t) - 1

		k1 = dy(t[j], y_rk4[j]) # this is just the euler's method (approximates the left)
		k2 = dy(t[j] + h/2, y_rk4[j] + (h/2)*k1) # the midpoint method
		k3 = dy(t[j] + h/2, y_rk4[j] + (h/2)*k2) # the midpoint method using the previous midpoint method
		k4 = dy(t[j] + h, y_rk4[j] + h*k3) # the full Euler's step (approximates the right)

		roc = (1/6)*(k1 + 2*k2 + 2*k3 + k4) # using the weighted average of the ks

		y_rk4[j + 1] = y_rk4[j] + h*roc # to euler's method using the weighted average
	end

	return y_rk4
end

##---------- System of ODEs Runge-Kutta ----------##

function rk4(F::Array{Function, 1}, t0::Float64, tn::Float64, Y0::Array{Float64, 1}, h::Float64, Y1::Function, Y2::Function, Y3::Function)
	
	t = collect(t0:h:tn)
	tlen = length(t)
	Y0len = length(Y0)
	# matrix that will store approximations for each function per time unit
	Y = zeros(Float64, tlen, Y0len)
	
	# vectors that will store the weighted rocs
	K1 = zeros(Float64, Y0len)
	K2 = zeros(Float64, Y0len)
	K3 = zeros(Float64, Y0len)
	K4 = zeros(Float64, Y0len)
	
	# at t = 0 (first row), all the Y values are the initial values (Y0s)
	for j = 1:Y0len Y[1,j] = Y0[j] end
	
	# same process as regular rk4
	for j = 1:tlen - 1
		
		for l = 1:Y0len K1[l] = F[l](t[j], Y[j,:]) end
		for l = 1:Y0len K2[l] = F[l](t[j] + h/2, Y[j,:] + K1.*h/2) end
		for l = 1:Y0len K3[l] = F[l](t[j] + h/2, Y[j,:] + K2.*h/2) end
		for l = 1:Y0len K4[l] = F[l](t[j] + h, Y[j,:] + K3.*h) end

		roc = (1/6).*(K1 + K2.*2 + K3.*2 + K4)

		Y[j + 1,:] = Y[j,:] + roc.*h
	end

	y1_tru = Y1.(t)
	y2_tru = Y2.(t)
	#y3_tru = Y3.(t)

	return ["t" "Y1" "Y2" "Y3" "Y1_tru" "Y2_tru";
	       t Y[:,1] Y[:,2] Y[:,3] y1_tru y2_tru]
end
