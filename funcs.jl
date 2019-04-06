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

##---------- Lagrange Polynomial Interpolation ----------##

# given some number of (x, y) points, interpolate through the points
# Takes: an array of x values and an array of y values
function L(x::Array{Float64, 1}, y::Array{Float64, 1}) # UNFINISHED! #
	
	n = length(x) # number of points given
	a = ones(n + 1) # degree of poly is n + 1 (at most)
	l = [[]] # uhhh.. lagrange things? not there yet

	# alternate-tively iterate through and slap on -1s
	# this is how the lagrange pieces are structured
	i = 2
	for j = 1:n + 1
		a[i] *= -1
		i += 2
		if(i > n + 1) break end
	end

	# meat of the code
	for i = 1:(n - 1)
		c = splice!(x, i) # take all but the current
		# the following are general for all polynomials
		# second term is the sum of all xs (except current)
		# the last term is the constant (product of all xs)
		a[2] = a[2]*reduce(+,x)
		a[end] = a[end]*reduce(*,x)
	
		# this only applies when the polynomail is greater
		# than degree 2. aka, given more than 3 points
		if n > 3
			tmp = []
			q = combos(x,tmp,1,2) # combinations of coefficients in pairs
			a[3] = a[3]*reduce(+,p) # standard 3rd term
			
			# loop from 4th term to last-1
			for j = 4:n
				for k = 1:(n - 1)
					for l = (n+1):length(q)
						a[j] = reduce(+,(x[k].*q[l:end]))
						
					end
				end
				print(a[j],"\n")
			end
		end
		append!([c],x) # add back the current
	end
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

##---------- Making life easier 1.0 ----------##

function Do(f::Function, P::Function, x::Array{Float64, 1}, xi::Function, a::Int64, b::Int64, g::Function, dg::Function)
	
	# finding abs err of lagrange poly
	print("f(x) = ", f(x), "\nP(x) = ", P(x),"\nAbs Err = ", abs(f(x) - P(x)))
	
	# finding err bounds of lagrange poly
	m = mindex(a, b, xi)
	M = maxdex(a, b, xi)
	print("\nMIN xi(x) = ", m[1]," at x = ", m[2], "\nMAX xi(x) = ", M[1], " at x = ",M[2])

	r = roots(Poly(dg))
	Gm = mindex(r, g)
	GM = maxdex(r, g)
	print("\nMIN g(r) = ", Gm[1], "\n\tfor g'(x) root: ", Gm[2], "\nMAX g(r) = ", GM[1], "\n\tfor g'(x) root: ", GM[2])
	
	print("\nErr Bounds: [", m[1]*Gm[1], ", ", M[1]*GM[1], "]")
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

##---------- Making Life Easier 2.0 ----------##

function Do(dy::Function, t0::Float64, tn::Float64, y0::Float64, h::Float64, tru::Function)
	
	# construct your time span
	t = collect(t0:h:tn)
	
	# allocate memory for your calculations
	y_eulr = zeros(Float64, length(t))
	y_mid = zeros(Float64, length(t))
	y_rk4 = zeros(Float64, length(t))
	
	# compute the actual solution and the approximations
	y_tru = tru.(t)
	y_eulr = eulr(dy, t, y0, y_eulr)
	y_mid = midpoint(dy, t, y0, y_mid)
	y_rk4 = rk4(dy, t, y0, y_rk4)
	
	# plot the actual and approximated solutions
	p0 = plot(t, [y_tru y_eulr y_mid y_rk4], title = "1a", label = ["y_tru" "y_eulr" "y_mid" "y_rk4"]);
	p1 = plot(t, y_tru); p2 = plot(t, y_eulr); p3 = plot(t, y_mid); p4 = plot(t, y_rk4);
	p5 = plot(p1, p2, p3, p4, label = ["y_tru" "y_eulr" "y_mid" "y_rk4"]);
	
	# compute the error between the approximations and the actual solution
	err_eulr = abs.(y_tru - y_eulr)
	err_mid = abs.(y_tru - y_mid)
	err_rk4 = abs.(y_tru - y_rk4)
	
	# plot the error
	p6 = plot(t, [err_eulr err_mid err_rk4], title = "3a", label = ["eulr" "mid" "rk4"]);
	
	savefig(p0,"fig1ai.png")
	savefig(p5,"fig1aii.png")
	savefig(p6,"fig3a.png")
	
	return ["t" "y_tru" "y_eulr" "err_eulr" "y_mid" "err_mid" "y_rk4" "err_rk4";
		t y_tru y_eulr err_eulr y_mid err_mid y_rk4 err_rk4]
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
