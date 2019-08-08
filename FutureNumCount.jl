using StatsBase

## PROBLEM: given an array of n numbers, [a, b, c, d, e, f] 
#           check into the future and find the amount of values greater than the current one
#           if c, e, f are all greater than a then the returning array will have a 3 in a's place
#           Example: given [2, 8, 1, 6, 7], the answer is [3, 0, 2, 1, 0]

# sample x
x = [2, 8, 1, 6, 7]
#x = [50, 13, 4, 76, 43, 12, 32, 21, 23, 87, 36, 9, 1, 5, 14, 11, 78, 64, 34]
#x = sample(collect(1:100000),100000,replace=false)

n = length(x) # just to avoid function calling length
y = zeros(n) # the returning answer array

function m(x) # finds max of x
	max = x[1]
	for j = 2:length(x) max = x[j] > max ? x[j] : max end
	return max
end

function s(nums, k, n) # calculate the sum of nums from k to length(nums)
	sum = 0
	for j = k:n sum += nums[j] end
	return sum
end

nums = zeros(m(x)) # so from 1 to n

y[n] = 0 # nothing is greater than the last value so 0
nums[x[n]] = 1 # already looked at last value, so just add 1 to it's index in nums (will only contain 1s)
N = length(nums) # length of nums (just to avoid function calling each time)

for j = collect(n-1:-1:1) # go backwards
	y[j] = s(nums, x[j]+1, N) # sum up all 1s from (the x[j] value + 1) till length of nums
	nums[x[j]] = 1 # add 1 to the x[j] index in nums
end

print(y) # look at results
