local F={}
function F.riemannC(f, n, a, b, init)
    expect(1, f, "function")
    expect(2, n, "number")
    expect(3, a, "number")
    expect(4, b, "number")
    expect(5, init, "number", "nil")
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

function F.normPDF(mu,sigma) return function(x) return math.exp(-0.5*((x-mu)/sigma)^2)/(sigma*math.sqrt(2*math.pi)) end end

--using normPDF, integrates from lower to upper with 10k subdivisions
function F.normalCDF(lower, upper, mean, stdev)
    expect(1, lower, "number")
    expect(2,upper,"number")
    expect(3,mean,"number")
    expect(4,stdev,"number")
    return F.riemannC(F.normPDF(x,mean,stdev),10000,lower,upper)
end
return F

print(F.normalCDF(-1, 1, 0, 1))