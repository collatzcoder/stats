--- A compilation and rework of the mmath file, incorporating optimized forms for statistics and newer math functions.
-- A rework, separated from the original cc:advanced math, so I can update it better
--by Mr James (sans.9536) of Atlas, dm me if you have questions 


local expect = require "cc.expect"
local expect = expect.expect
local F = {}


--{TABLE OPS}

--linear search for whatever value you want
function F.find(tbl, value)
    expect(1, tbl, "table")
    for i=1, #tbl do
        if tbl[i] == value then
            return i
        end
    end
    return nil
end

--make n numbers between a and b
function F.testdata(n, a, b)
    expect(1, n, "number")
    expect(2, a, "number")
    expect(3, b, "number")
    if n < 1 then expect(1, n, "number greater than or equal to one")
    local result = {}
    a = a or 0
    b = b or 1
    for i = 1, n do
        result[i] = a + math.random() * (b - a)
    end
    return result
end
--add each value in values to the resultant table n times
function F.weightedTable(values, weights)
    expect(1, values, "table")
    expect(2, weights, "table")
    if #values ~= #weights then
        error("The values and weights tables must have the same number of values!")
    end

    local result = {}
    for i,v in pairs(values) do
        local b = 1
        while weights[i] >= b do
            table.insert(result,v)
            b = b + 1
        end
    end
    return result
end

--take each value in the table and insert it into a random position in the table, effectively scrambling it
--not perfect but if you ask me it's good enough
function F.scramble(t)
    expect(1, t, "table")
    for index, value in pairs(t) do
        local a = value
        table.remove(t, index)
        table.insert(t, math.random(1, #t + 1), a)
    end
    return t
end

--{FUNCTIONS}
--WIP



--{CALCULUS}

--[DERIVATIVES]
--computes average or instantaneous rate of change(returns a value). shrink h to approach the IRC, defaults to 1e-10
--tested to be accurate to 1e-5
function F.AIRC(f, x, h)
    expect(1, f, "function")
    expect(2, x, "number")
    expect(3, h, "number", "nil")

    h = h or 1e-10
    return (f(x + h) - f(x - h)) / (2 * h)
end
--returns a function that represents the derivative of f(x), will try to make something that builds an actual function that can be printed to a string
function F.ddx(f)
    expect(1, f, "function")
    return function(x)
        return F.AIRC(f, x)
    end
end


--[INTEGRALS]
--trapezoidal riemann sum with n subdivisions from a to b
function F.riemannS(f, n, a, b, init)
    expect(1, f, "function")
    expect(2, n, "number")
    expect(3, a, "number")
    expect(4, b, "number")
    expect(5, init, "number", "nil")
    if n <= 0 then
        error("Number of Subdivisions must be greater than zero!")
    end
    if a >= b then
        error("The lower bound must be less than the upper bound!")
    end

    local width = (b - a) / n
    local sum = init or 0
    for i = 0, n - 1 do
        local x1 = a + i * width
        local x2 = a + (i + 1) * width
        sum = sum + (f(x1) + f(x2)) * width / 2
    end
    return sum
    
end

--use simpsons rule to integrate with n subdivisions from a to b
--MAY NOT WORK I THINK ITS BROKEN
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
--note: you cannot have an integral at a point. I cannot program something to return the function that is the integral (yet)

--{SOLVERS}
--uses newton raphson method to find the root, can fail if the derivative is zero or near zero, or if it fails to converge within 500 iterations
function F.solveRoot(f, x0, tol)
    expect(1, f, "function")
    expect(2, x0, "number")
    expect(3, tol, "number")

    local x = x0
    for i = 1, 500 do
        local dfx = F.ARC(f, x)
        if math.abs(dfx) < 1e-12 then
            return nil
        end
        local x_new = x - f(x) / dfx
        if math.abs(x_new - x) < tol then
            return x_new
        end
        x = x_new
    end
    return nil
end
--solve system of equations by finding solutions of f-g
function F.solveSysEq(f, g, min, max, steps, tol, close_thresh)
    expect(1, f, "function")
    expect(2, g, "function")
    expect(3, min, "number")
    expect(4, max, "number")
    expect(5, steps, "number")
    expect(6, tol, "number", "nil")
    expect(7, close_thresh, "number", "nil")
    if min >= max then
        error("Minimum bound needs to be lesser than maximum bound!")
    end
    if steps < 1 then
        error("Number of Steps must be greater than or equal to one!")
    end
    local function h(f,g) 
        return function(x) return (f(x)-g(x)) end
    end
    local solutions = {}
    for i = 0, steps do
        local root = F.solveRoot(h(f,g), min + (max - min) * i / steps, tol)
        if root then
            local duplicate = false
            for _, sol in ipairs(solutions) do
                if math.abs(sol - root) < close_thresh then
                    duplicate = true
                    break
                end
            end
            if not duplicate and root >= min and root <= max then
                table.insert(solutions, root)
            end
        else
            local x = min + (max - min) * i / steps
            for j = 1,100 do
                if h(f,g)(x) < close_thresh then
                    local duplicate = false
                    for _, sol in ipairs(solutions) do
                        if math.abs(sol - x) < close_thresh then
                            duplicate = true
                            break
                        end
                    end
                    if not duplicate then
                        table.insert(solutions, x)
                    end
                    break
                else
                    x = x + (max - min) / (steps * 100)
                end
            end
        end
    end
    return solutions
end

--{STATISTICS}
--
-- @section basic_descriptive_statistics

--- Computes the sum of numeric values in the dataset
--
-- @tparam number[] data sequential table of numbers
-- @treturn number sum of values (0 if empty)
function F.sum(data)
    expect(1, data, "table")
    local s = 0
    for _, v in pairs(data) do
        if type(v) ~= "number" then expect(1, data, "sequential table of numbers") end
        s = s + v
    end
    return s
end

--- Computes the arithmetic mean (average) of the dataset
--
-- @tparam number[] data
-- @treturn number|nil mean or nil if data empty
function F.mean(data)
    local n = #data
    return n > 0 and F.sum(data) / n or nil
end

--- Computes the minimum value of the dataset
--
-- @tparam number[] data
-- @treturn number|nil minimum value or nil if data empty
function F.min(data)
    expect(1, data, "table")
    table.sort(data)
    return data[1]
end

--- Computes the maximum value of the dataset
--
-- @tparam number[] data
-- @treturn number|nil maximum value or nil if data empty
function F.max(data)
    expect(1, data, "table")
    table.sort(data)
    return data[#data]
end

--- Computes population variance by default.
--
-- Set `isSample` to true for sample variance.
-- @tparam number[] data
-- @tparam[opt=false] boolean isSample true to compute sample variance (divide by n-1)
-- @treturn number|nil variance or nil if not defined
function F.variance(data, isSample)
    expect(1, data, "table")
    if is_sample ~= nil then expect(2, isSample, "boolean") end
    if #data == 0 then return nil end
    local sumSq = 0
    for _, v in pairs(data) do
        sumSq = sumSq + (v - F.mean(data)) ^ 2
    end
    local denom = (isSample and #data - 1 or #data)
    if denom <= 0 then return nil end
    return sumSq / denom
end

--- Computes the standard deviation of the dataset using sample variance
--
-- @tparam number[] data
-- @tparam[opt=false] boolean isSample true to compute sample stdev
-- @treturn number|nil standard deviation or nil if not defined
function F.stdev(data, isSample)
    local v = F.variance(data, isSample)
    return v and math.sqrt(v) or nil
end

--- Number of elements in the dataset
--  If you don't know how this works i hate you
--
-- @tparam table data
-- @treturn number size (n)
function F.size(data)
    expect(1, data, "table")
    return #data
end

--- Computes the range of the dataset(max-min). Keep in mind that this is susceptible to skewing by outliers.
--
-- @tparam number[] data
-- @treturn number|nil range or nil if undefined
function F.range(data)
    return F.min(data) and F.max(data) and F.max(data) - F.min(data) or nil
end

--- Computes the median (50th percentile) by sorting then popping one off of each side
--
-- If even length, returns average of the two middle values
-- @tparam number[] data
-- @treturn number|nil median or nil if empty
function F.median(data)
    expect(1, data, "table")
    table.sort(data)
    repeat
        table.remove(data,1)
        table.remove(data)
    until #data<=2
    return F.sum(data)/#data
end

--- Computes the mode and returns a sorted array of the most frequent value(s)
--
-- @tparam table data
-- @treturn table list of mode values (empty table if data empty)
function F.mode(data)
    expect(1, data, "table")
    local counts = {}
    local maxc = 0

    for _, v in pairs(data) do
        counts[v] = (counts[v] or 0) + 1
        if counts[v] > maxc then
            maxc = counts[v]
        end
    end
    if maxc == 0 then
        return {}
    end

    local modes = {}
    for v, c in pairs(counts) do
        if c == maxc then
            table.insert(modes, v)
        end
    end
    table.sort(modes)
    return modes
end

--- Computes the standard error of the mean (commonly called SE) indicating precision
--
-- @tparam number[] data
-- @treturn number|nil standard error or nil if undefined
function F.sex(data)
    return #data > 0 and F.stdev(data, true) / math.sqrt(#data) or nil
end
--- Computes the asymmetry (skewness) of data around its mean
--
-- Requires at least 3 observations
-- @tparam number[] data
-- @treturn number|nil skewness or nil if undefined
function F.skewness(data)
    expect(1, data, "table")
    local n = #data
    if n < 3 then return nil end
    local mu, sd = F.mean(data), F.stdev(data, true)
    if not sd or sd == 0 then return nil end
    local m3 = 0
    for _, v in pairs(data) do
        m3 = m3 + ((v - mu) / sd) ^ 3
    end
    return n / ((n - 1) * (n - 2)) * m3
end

--- Computes the messure of skew (kurtosis) using an excess kurtosis formula variant
--
-- Requires at least 4 observations and indicates how much peak or tail a distribution would have
-- @tparam number[] data
-- @treturn number|nil kurtosis or nil if undefined
function F.kurtosis(data)
    expect(1, data, "table")
    local n = #data
    if n < 4 then return nil end
    local mu, sd = F.mean(data), F.stdev(data, true)
    if not sd or sd == 0 then return nil end
    local m4 = 0
    for _, v in pairs(data) do
        m4 = m4 + ((v - mu) / sd) ^ 4
    end
    return (n * (n + 1) / ((n - 1) * (n - 2) * (n - 3))) * m4
        - (3 * (n - 1) ^ 2) / ((n - 2) * (n - 3))
end

--- Alternative Means
--
-- @section alternative_means

--- Computes the geometric mean
--
-- Only defined for strictly positive values
--
-- Useful for multiplicative rates such as combining percentage growth rates
-- @tparam number[] data
-- @treturn number|nil geometric mean or nil if any value <= 0 or data empty
function F.geometricMean(data)
    expect(1, data, "table")
    local prod, n = 1, #data
    if n == 0 then return nil end
    for _, v in pairs(data) do
        if type(v) ~= "number" or math.abs(v) ~= v then expect(1, data, "sequential table of positive numbers") end
        if v <= 0 then return nil end
        prod = prod * v
    end
    return prod ^ (1 / n)
end

--- Computes the harmonic mean
--
-- Best used for average rates and ratios, or inverse proportionalities
-- @tparam number[] data
-- @treturn number|nil harmonic mean or nil if any value == 0 or data empty
function F.harmonicMean(data)
    expect(1, data, "table")
    local sumr, n = 0, #data
    if n == 0 then return nil end
    for _, v in pairs(data) do
        if type(v) ~= "number" or v == 0 or math.abs(v) ~= v then expect(1, data, "sequential table of positive non-zero numbers") end
        if v == 0 then return nil end
        sumr = sumr + 1 / v
    end
    return n / sumr
end

--- Computes the trimmed mean
--
-- Removes `trim` fraction from each tail and averages the rest
--
-- Much more reliable measure of central tendency, clamps down outliers by removing a percentage of extremes
-- @tparam number[] data
-- @tparam[opt=0.05] number trim fraction to trim from each tail
-- @treturn number|nil trimmed mean or nil if undefined
function F.trimmedMean(data, trim)
    if trim ~= nil then
        expect(2, trim, "number")
        if math.abs(trim) ~= trim then
            expect(2, trim, "non-negative number")
        end
    end
    trim = trim or 0.05
    local t = { table.unpack(data) }
    table.sort(t)
    local n = #t
    if n == 0 then return nil end
    local tn = math.floor(trim * n)
    local s = 0
    for i = tn + 1, n - tn do
        if type(t[i]) ~= "number" then expect(1, data, "sequential table of numbers") end
        s = s + t[i]
    end
    local denom = n - 2 * tn
    if denom <= 0 then return nil end
    return s / denom
end

--- Computes the winsorized mean
--
-- Clamps extreme values to given quantiles before averaging
--
-- Much more reliable measure of central tendency
-- @tparam number[] data
-- @tparam[opt=0.05] number alpha fraction to winsorize in each tail
-- @treturn number|nil winsorized mean or nil if empty
function F.winsorMean(data, alpha)
    if alpha ~= nil then
        expect(2, alpha, "number")
        if math.abs(alpha) ~= alpha then
            expect(2, alpha, "non-negative number")
        end
        if alpha > 0.5 then execpt(2, alpha, "number less than or equal to 0.5")
    end
    alpha = alpha or 0.05
    local t = { table.unpack(data) }
    table.sort(t)
    local n = #t
    if n == 0 then return nil end
    if alpha == 0 then return F.mean(t) end

    local lower = math.floor(alpha * n) + 1
    local upper = math.ceil((1 - alpha) * n)
    local winsorT = { table.unpack(t) }

    local lowerBound = t[lower]
    local upperBound = t[upper]
    for i = 1, lower - 1 do winsorT[i] = lowerBound end
    for i = upper + 1, n do winsorT[i] = upperBound end

    return F.mean(winsorT)
end

--- Computes the weighted mean
--
-- @tparam number[] data
-- @tparam number[] weights same length as data
-- @treturn number|nil weighted mean or nil if lengths mismatch or zero total weight
function F.weightedMean(data, weights)
    expect(1, data, "table")
    expect(2, weights, "table")
    local n = #data
    if n == 0 or #weights ~= n then return nil end
    local num, den = 0, 0
    for i = 1, n do
        if type(data[i]) ~= "number" then expect(1, data, "sequential table of numbers") end
        if type(weights[i]) ~= "number" then expect(2, weights, "sequential table of numbers") end
        num = num + data[i] * weights[i]
        den = den + weights[i]
    end
    return den > 0 and num / den or nil
end

--- Robust Statistics
--
-- @section robust_statistics

--- Computes the median absolute deviation (MAD)
--
-- Much more reliable measure of variability when outliers are present, from the deviation
-- @tparam number[] data
-- @treturn number|nil MAD or nil if empty
function F.madMedian(data)
    if #data == 0 then return nil end
    local med = F.median(data)
    local devs = {}
    for _, v in pairs(data) do
        if type(v) ~= "number" then expect(1, data, "sequential table of numbers") end
        table.insert(devs, math.abs(v - med))
    end
    return F.median(devs)
end

--- Inequality Measures
--
-- @section inequality_measures

--- Computes the gini coefficient for inequality (0..1)
--
-- Returns 0 for datasets with fewer than 2 or zero mean
--
-- Represents how unequally distributed a dataset it, best explained by income but works for other sets too
-- @tparam number[] data
-- @treturn number gini coefficient (0..1)
function F.gini(data)
    expect(1, data, "table")
    local n = #data
    if n < 2 then return 0 end
    local t = { table.unpack(data) }
    table.sort(t)
    local mean = F.mean(t)
    if not mean or mean == 0 then return 0 end

    local cumsum = 0
    for i = 1, n do
        if type(t[i]) ~= "number" then expect(1, data, "sequential table of numbers") end
        cumsum = cumsum + t[i]
    end

    local num = 0
    for i = 1, n do
        num = num + (2 * i - n - 1) * t[i]
    end
    return num / (n * cumsum)
end

--- Bivariate Statistics
--
-- @section bivariate_statistics

--- Computes the sample covariance (uses n-1 denominator)
--
-- Represents how much the two datasets will vary with each other
-- @tparam number[] x
-- @tparam number[] y
-- @treturn number|nil covariance or nil if lengths mismatch or n<=1
function F.covariance(x, y)
    expect(1, x, "table")
    expect(2, y, "table")
    local n = #x
    if n == 0 or #y ~= n then return nil end
    local mx, my = F.mean(x), F.mean(y)
    local cov = 0
    for i = 1, n do
        if type(x[i]) ~= "number" then expect(1, x, "sequential table of numbers") end
        if type(y[i]) ~= "number" then expect(2, y, "sequential table of numbers") end
        cov = cov + (x[i] - mx) * (y[i] - my)
    end
    if n <= 1 then return nil end
    return cov / (n - 1)
end

--- Computes the pearson correlation coefficient
--
-- Returns 0 if undefined
--
-- Represents the strength and relationship of the covariances, much more interpretable
-- @tparam number[] x
-- @tparam number[] y
-- @treturn number correlation coefficient in [-1,1] or 0 if undefined
function F.correlation(x, y)
    local cov = F.covariance(x, y)
    local sx, sy = F.stdev(x, true), F.stdev(y, true)
    if not cov or not sx or not sy or sx == 0 or sy == 0 then
        return 0
    end
    return cov / (sx * sy)
end

--- Quantiles and Outliers
--
-- @section quantiles_and_qutliers

--- Computes the percentile interpolation (linear interpolation between order statistics)
--
-- Finds the percentile of a value p given a dataset, sorts and places in between indices(or on one)
--
-- p in [0,1]
-- @tparam number[] data
-- @tparam number p percentile proportion (0..1)
-- @treturn number|nil percentile value or nil if data empty
function F.percentile(data, p)
    expect(2, p, "number")
    if p < 0 then expect(2, p, "number greater than or equal to zero") end
    if p > 1 then expect(2, p, "number less than or equal to one") end
    if p == 0 then return F.min(data) end
    if p == 1 then return F.max(data) end
    local t = { table.unpack(data) }
    table.sort(t)
    local n = #t
    if n == 0 then return nil end
    local rank = p * (n - 1) + 1
    local f, c = math.floor(rank), math.ceil(rank)
    if type(t[f]) ~= "number" then expect(1, data, "sequential table of numbers") end
    if f == c then
        return t[f]
    else
        if type(t[c]) ~= "number" then expect(1, data, "sequential table of numbers") end
        return t[f] + (rank - f) * (t[c] - t[f])
    end
end

--- Computes the first and third quartiles as table (`{ Q1 = ..., Q3 = ... }`)
--
-- Quartiles are found by cutting the dataset in half using the median, and finding the median of those sets
-- @tparam number[] data
-- @treturn table quartiles
function F.quartiles(data)
    return {
        Q1 = F.percentile(data, 0.25),
        Q3 = F.percentile(data, 0.75)
    }
end

--- Computes the interquartile range Q3 - Q1
--
-- The distance between first and third quartiles
-- @tparam number[] data
-- @treturn number|nil IQR or nil if undefined
function F.iqr(data)
    local q = F.quartiles(data)
    if not q.Q1 or not q.Q3 then return nil end
    return q.Q3 - q.Q1
end

--- Computes a list of all outliers in the dataset
--
-- Outliers are considered as such if they are more than 150% of the IQR away from either the first or third quartiles
-- @tparam number[] data
-- @treturn number[] list of outlier values (empty if none)
function F.outliers(data)
    local q = F.quartiles(data)
    local iqr = F.iqr(data)
    if not iqr then return {} end
    local low, high = q.Q1 - 1.5 * iqr, q.Q3 + 1.5 * iqr
    local outs = {}
    for _, v in pairs(data) do
        if v < low or v > high then
            table.insert(outs, v)
        end
    end
    return outs
end

--- Utility Functions


--- Linear Regression
--
-- @section linear_regression

--- Computes a simple least-squares linear regression (y ~ a + b x)
--
-- Minimizes squares to find the line that is the least distance squared from all data points
-- @tparam number[] x independent variable values
-- @tparam number[] y dependent variable values
-- @treturn table|nil `{ slope = number, intercept = number }` or nil if undefined
function F.linReg(x, y)
    expect(1, x, "table")
    expect(2, y, "table")
    local n = #x
    if n == 0 or #y ~= n then return nil end

    local mx, my = F.mean(x), F.mean(y)
    local num, den = 0, 0
    for i = 1, n do
        if type(x[i]) ~= "number" then expect(1, x, "sequential table of numbers") end
        if type(y[i]) ~= "number" then expect(2, y, "sequential table of numbers") end
        local dx = x[i] - mx
        num = num + dx * (y[i] - my)
        den = den + dx * dx
    end
    if den == 0 then return nil end

    local slope = num / den
    local intercept = my - slope * mx
    return { slope = slope, intercept = intercept }
end

--- Predicts the y value for a given linear regression model at a given position
--
-- @tparam table model `{ slope, intercept }`
-- @tparam number xval value of x
-- @treturn number|nil predicted y or nil if model missing
function F.linRegPred(model, xval)
    expect(1, model, "table")
    expect(2, xval, "number")
    expect(1, model.slope, "linear regression model")
    expect(1, model.intercept, "linear regression model")
    if not model then return nil end
    return model.slope * xval + model.intercept
end

--- Computes the coefficient of determination R^2 for linear model
--
-- Represents how much variation in variable y can be attributed to variable x
-- @tparam number[] x
-- @tparam number[] y
-- @tparam table model linear model returned by linReg
-- @treturn number|nil R^2 (1 if total variance is zero)
function F.r2(x, y, model)
    expect(1, x, "table")
    expect(2, y, "table")
    expect(3, model, "table")
    expect(3, model.slope, "linear regression model")
    expect(3, model.intercept, "linear regression model")
    local n = #y
    if n == 0 or not model then return nil end

    local ymean = F.mean(y)
    local ssres, sstot = 0, 0
    for i = 1, n do
        if type(x[i]) ~= "number" then expect(1, x, "sequential table of numbers") end
        if type(y[i]) ~= "number" then expect(2, y, "sequential table of numbers") end
        local ypred = F.linRegPred(model, x[i])
        ssres = ssres + (y[i] - ypred) ^ 2
        sstot = sstot + (y[i] - ymean) ^ 2
    end
    if sstot == 0 then
        return 1
    else
        return 1 - ssres / sstot
    end
end



--returns the function for the normal curve
function F.normPDF(mu,sigma) return function(x) return math.exp(-0.5*((x-mu)/sigma)^2)/(sigma*math.sqrt(2*math.pi)) end end

--using normPDF, integrates from lower to upper with 10k subdivisions
function F.normCDF(lower, upper, mean, stdev)
    expect(1, lower, "number")
    expect(2,upper,"number")
    expect(3,mean,"number")
    expect(4,stdev,"number")
    return F.riemannC(F.normPDF(x,mean,stdev),10000,lower,upper)
end

--for now, normal gamma function, only takes positive integers
function F.gammaN(x) return math.factorial(x-1) end

function F.tPDF(mu,sigma,df)



return F