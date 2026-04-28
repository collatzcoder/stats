local function riemannC(f, n, a, b, init)
    if n <= 2 then
        error("Number of Subdivisions must be greater than two!")
    end
    if a >= b then
        error("The lower bound must be less than the upper bound!")
    end

    n=math.floor(n/2)*2 
    local width = (b - a) / n
    local sum = init or 0
    for i = 0, n / 2 - 1 do
        local x0 = a + i * 2 * width
        local x1 = a + (i * 2 + 1) * width
        local x2 = a + (i * 2 + 2) * width
        sum = sum + (f(x0) + 4 * f(x1) + f(x2)) * width / 3
    end
    return sum
    
end

local function integral(f)
    
    return function(x)
        return riemannC(f,1e5,x,x+1e-2)
    end
end
print(integral(function(x) return x^2 end)(2)) -- should be close to 8/3