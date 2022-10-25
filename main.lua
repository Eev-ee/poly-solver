

local complex
local complex_operation = {}
complex_operation.__mul = function(lhs, rhs)
	lhs = type(lhs) == "table" and lhs or type(lhs) == "number" and complex(lhs, 0) or error("MEOW!!")
	rhs = type(rhs) == "table" and rhs or type(rhs) == "number" and complex(rhs, 0) or error("MEOW!!")

	local a, b, c, d = lhs.re, lhs.im, rhs.re, rhs.im
	return complex(a*c - b*d, b*c + a*d, true)
end
complex_operation.__div = function(lhs, rhs)
	lhs = type(lhs) == "table" and lhs or type(lhs) == "number" and complex(lhs, 0) or error("MEOW!!")
	rhs = type(rhs) == "table" and rhs or type(rhs) == "number" and complex(rhs, 0) or error("MEOW!!")

	local a, b, c, d = lhs.re, lhs.im, rhs.re, rhs.im
	local m = rhs:abs2()
	return complex((a*c + b*d)/m, (b*c - a*d)/m, true)
end
complex_operation.__add = function(lhs, rhs)
	lhs = type(lhs) == "table" and lhs or type(lhs) == "number" and complex(lhs, 0) or error("MEOW!!")
	rhs = type(rhs) == "table" and rhs or type(rhs) == "number" and complex(rhs, 0) or error("MEOW!!")

	local a, b, c, d = lhs.re, lhs.im, rhs.re, rhs.im
	return complex(a + c, b + d, true)
end
complex_operation.__sub = function(lhs, rhs)
	lhs = type(lhs) == "table" and lhs or type(lhs) == "number" and complex(lhs, 0) or error("MEOW!!")
	rhs = type(rhs) == "table" and rhs or type(rhs) == "number" and complex(rhs, 0) or error("MEOW!!")

	local a, b, c, d = lhs.re, lhs.im, rhs.re, rhs.im
	return complex(a - c, b - d, true)
end
complex_operation.__pow = function(lhs, rhs)
	lhs = type(lhs) == "table" and lhs or type(lhs) == "number" and complex(lhs, 0) or error("MEOW!!")
	rhs = type(rhs) == "table" and rhs or type(rhs) == "number" and complex(rhs, 0) or error("MEOW!!")

	local a, b, c, d = lhs.re, lhs.im, rhs.re, rhs.im
	local arg = lhs:arg()
	local k = lhs:abs2()
	local lk = math.log(k)
	local m = k^(c/2)*math.exp(-d*arg)
	local t = c*arg + d/2*lk
	return complex(math.cos(t)*m, math.sin(t)*m, true)
end
complex_operation.__tostring = function(self)
	local re, im = self.re, self.im
	--re, im = re - re%0.0001, im - im%0.0001
	return string.format("%.20f", re) .. (im < 0 and " - " .. string.format("%.20f", -im) or " + " .. string.format("%.20f", im)) .. "I"
end

function complex( re, im, filter )
	if math.abs(im) <= 1e-10 and filter then
		return re
	end
	return setmetatable({
		re = re,
		im = im,
		abs = function(self)
			return (re*re + im*im)^0.5
		end,
		abs2 = function(self)
			return re*re + im*im
		end,
		arg = function(self)
			return math.atan2(im, re)
		end
	}, complex_operation)
end

print(complex(1, 2) + 2)


local function p( coeff_array, x )
	local result = 0
	for i = 1, #coeff_array do
		result += coeff_array[i]*x^(i - 1)
	end
	return result
end
local function p_d( coeff_array, x )
	local n = #coeff_array
	if n == 1 then
		return 0
	end
	local result = 0
	for i = 1, n - 1 do
		result += i*coeff_array[i + 1]*x^(i - 1)
	end
	return result
end
local function p_dd( coeff_array, x )
	local n = #coeff_array
	local result = 0
	for i = 3, n do
		result += (i - 1)*(i - 2)*x^(i - 3)*coeff_array[i]
	end
	return result
end
local function p_ddd( coeff_array, x )
	local n = #coeff_array
	local result = 0
	for i = 4, n do
		result += (i - 1)*(i - 2)*(i - 3)*x^(i - 4)*coeff_array[i]
	end
	return result
end

local function div_diff(coeff_array, xs)
	local n = #xs
	if n == 2 then
		local a, b = xs[1], xs[2]
		return (p(coeff_array, b) - p(coeff_array, a))/(b - a)
	else
		local a = {}
		local b = {}
		for i = 2, n do
			a[i - 1] = xs[i]
		end
		for i = 1, n - 1 do
			b[i] = xs[i]
		end
		return (div_diff(coeff_array, a) - div_diff(coeff_array, b))/(xs[n] - xs[1])
	end
end

local function symbolic_p_d( coeff_array )
	local result = {}
	for i = 1, #coeff_array - 1 do
		result[i] = i*coeff_array[i + 1]
	end
	return result
end

local function sqrt( x )
	if type(x) == "number" then
		if x < 0 then
			return complex(0, (-x)^0.5)
		end
		return x^0.5
	end
	local re, im = x.re, x.im
	local m = x:abs()
	return complex(((m + re)/2)^0.5, ((m - re)/2)^0.5*(im >= 0 and 1 or -1))
end

local function laguerre( coeff_array, initial_guess, max_iter, precision )
	precision = precision or 1e-5
	max_iter = max_iter or 1/0

	local guess = initial_guess or 0
	local n = #coeff_array - 1
	local iter = 0
	while true do
		iter += 1
		if iter >= max_iter then break end
		local px = p(coeff_array, guess)
		if type(px) == "table" and px.re then
			if math.abs(px.re) <= precision and math.abs(px.im) <= precision then break end
		elseif math.abs(px) <= precision then
			break
		end

		local G = p_d(coeff_array, guess)/px
		local H = G*G - p_dd(coeff_array, guess)/px
		local a = n/(G + sqrt((n - 1)*(n*H - G*G)))
		guess -= a
	end
	return guess 
end

local function householder_third( coeff_array, initial_guess, max_iter, precision )
	precision = precision or 1e-5
	max_iter = max_iter or 1/0

	local guess = initial_guess or complex(1, 1, false)
	for i = 1, max_iter do
		local fx = p(coeff_array, guess)
		if type(fx) == "table" and fx.re then
			if math.abs(fx.re) <= precision and math.abs(fx.im) <= precision then break end
		elseif math.abs(fx) <= precision then
			break
		end
		local fpx = p_d(coeff_array, guess)
		local fppx = p_dd(coeff_array, guess)
		local fpppx = p_ddd(coeff_array, guess)
		guess = guess - (6*fx*fpx*fpx - 3*fx*fx*fppx)/(6*fpx*fpx*fpx - 6*fx*fpx*fppx + fx*fx*fpppx)
	end
	return guess
end

local function muller( coeff_array, initial_guess, max_iter, precision )
	precision = precision or 1e-5
	max_iter = max_iter or 1/0

	local x0, x1, x2 = initial_guess - 1, initial_guess, initial_guess + 1
	for i = 1, max_iter do
		local w = div_diff(coeff_array, {x2, x1}) + div_diff(coeff_array, {x2, x0}) + div_diff(coeff_array, {x2, x1})
		local s_delta = sqrt(w*w - 4*p(coeff_array, x2)*div_diff(coeff_array, {x2, x1, x0}))
		local x3 = x2 - 2*p(coeff_array, x2)/(w + s_delta)
		if type(x3) == "table" and x3.re then
			if math.abs(x3.re) <= precision and math.abs(x3.im) <= precision then break end
		elseif math.abs(x3) <= precision then
			break
		end
		x0, x1, x2 = x1, x2, x3
	end
	return x2
end

local function solve_poly( coeff_array, initial_guess, max_iter, precision, result )
	if #coeff_array > 1 then
		--local root = laguerre(coeff_array, initial_guess, max_iter, precision)
		--local root = muller(coeff_array, initial_guess, max_iter, precision)
		local root = householder_third(coeff_array, initial_guess, max_iter, precision)
		local new_coef = {coeff_array[#coeff_array]}
		for i = #coeff_array - 1, 2, -1 do
			table.insert(new_coef, 1, coeff_array[i] + new_coef[1]*root)
		end
		table.insert(result, root)
		solve_poly(new_coef, root, max_iter, precision, result)
	else
		return
	end
end

local function print_poly( coeff )
	local result = {}
	for i = 1, #coeff do
		local l = coeff[i]
		if l == 0 then continue end
		if i > 1 then
			table.insert(result, l < 0 and "-" or "+")
			table.insert(result, math.abs(l) .. "x^" .. (i - 1))
		else
			table.insert(result, l)
		end
	end
	return table.concat(result, " ")
end

local function kahan_babuska(input)
	local sum = 0.0
	local cs = 0.0
	local ccs = 0.0
	local c = 0.0
	local cc = 0.0

	for i = 1, #input do
		local t = sum + input[i]
		if math.abs(sum) >= math.abs(input[i]) then
			c = (sum - t) + input[i]
		else
			c = (input[i] - t) + sum
		end
		sum =  t
		t = cs + c
		if math.abs(cs) >= math.abs(c) then
			cc = (cs - t) + c
		else
			cc = (c - t) + cs
		end
		cs = t
		ccs = ccs + cc
	end

	return sum + cs + ccs
end



local real = os.clock()
local result = {}
local coeff = {0, -36, 0, 49, 0, -14, 0, 1}
print(print_poly(coeff))
solve_poly(coeff, complex(1, 1, false), 500, 1e-10, result)
warn(os.clock() - real)
local errors = {}
for i = 1, #result do
	errors[i] = p(coeff, result[i])
	errors[i] = type(errors[i]) == "table" and errors[i]:abs() or math.abs(errors[i])
end

local error_avg = kahan_babuska(errors)/#errors
print(error_avg)
warn(kahan_babuska({0.1, 0.1, 0.1}))
warn("roots : ")
table.foreach(result, warn)
