function intPartition(n::Int; k::Int = 1)
	
	# Base Case: the only partion we actually know the answer to
	if n == 1 return [[1]] end

	# Otherwise:
	retval = [[n]] # for sure will include n itself
	
	# goes up to (n-k) this is to avoid duplicates
	# past this threshold, the partitions are repeated in reverse orders
	while k <= n-k # can also do n/2 works out to be the same thing
		
		# find all the partitions that would include k, where n = (n-k)
		# to find all the smaller sub partitions
		sublist = intPartition(n-k, k=k)

		# now that we have the subpartitions for THAT k, we want to tag
		# k onto that sublist to complete that parition for n
		map(sublistLists -> prepend!(sublistLists, k), sublist)

		# now append the found paritions to the retval
		append!(retval, sublist)

		k += 1
	end

	return retval
end
