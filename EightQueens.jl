# The recursive function that actually finds and places 
# the queens in a possible spot
function EQ(A, row)
	
	# Base Case: if you've reached a row that doesn't exist,
	# you've completely placed all possible queens
	# returning 1 adds to count
	if row > size(A)[1] return 1
	end
	
	# Initializing count at 0
	# (for every row, it'll restart to count,
	# but at every return, the final value of count is added
	# to the first call ever made)
	count = 0
	
	# The meat of the problem: you traverse every column
	# till you've reached the end, you check if you can place
	# a queen, if you can, then do so, label it as 1 so you are
	# aware of where you cannot place for the following rows
	for col = 1:size(A)[1]
		
		if chk(A, row, col) == true
			A[row, col] = 1
			count += EQ(A, row + 1) # once this comes back, you add
						# it to count, row now goes back to
						# what it was in the previous call
			A[row, col] = 0 	# you want to clear the placement of
						# the queen on that row, so you can find 
						# all the other queen placements in that row
		end
	end
	return count
end

# This is the checker that checks (how surprising o_o)
function chk(A, row, col)
	
	i = row # i & j is just
	j = col # for reuseability
	while i > 0 && j > 0 # i AND j decrease (basically: checking ALL of the major diagonal)
		
		if A[i, j] == 1 # Just don't place if there's something there already. Just don't do it lol
			return false
		end
		i -= 1
		j -= 1
	end

	i = row # The actual act
	j = col # of reusing :]
	
	while i > 0 # only decreasing i (checking all of the above, seriously only up)
		
		if A[i, j] == 1
			return false
		end
		i -= 1
	end

	i = row
	j = col
	
	while i > 0 && j <= size(A)[1] # decreasing i but INCREASING j (this checks the minor diagonal)
			
		if A[i, j] == 1
			return false
		end
		i -= 1
		j += 1
	end
	return true # if at the end of it all, nothing was in the way, you can place a queen :D
end

n = 10
A = zeros(Int64,n,n)
@time count = EQ(A, 1)
println(count, " solutions for a ", n, "x",n," matrix")
