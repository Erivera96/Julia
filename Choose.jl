function choose(xs::Vector, r::Int)

	# Base Case: if r = 1, want to return a list of list
	# because it's the only thing we for sure know how to do
	if r == 1 return [[x] for x in xs] end
	
	# else get all the combinations components 
	retval = []
	
	# don't want to go to end because doesn't make sense to
	# append smaller than r
	for i = 1:length(xs)-r+1
		
		# get smaller solutions for each iteration
		sublist = choose(xs[i+1:end], r-1)

		# put it together with CURRENT iteration
		map(sublistLists -> prepend!(sublistLists, xs[i]), sublist)
		
		# glue solutions together
		append!(retval, sublist)
	end
	
	# return list with all combinations
	return retval
end
