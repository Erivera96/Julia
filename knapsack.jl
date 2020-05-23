function knpsk(Capacity::Int64, numItems::Int64, items::Vector)
	# Note: Capacity is the maximum weight the knapsack can hold
	# numItems is the number of items available to choose from and 
	# items is given in (value, weight)

	tab = zeros(Int64, numItems+1, Capacity+1) # will store previous item's best

	for item = 2:numItems+1 # tab is padded with 0s, so items[1][1] = tab[2][2]
		for weight = 1:Capacity
			if items[item-1][2] <= weight # if the current item fits in the current weight limit then:
				# case 1: check if putting the item plus the best of the remaining
				# weight is better than the best of the previous same weight
				tab[item, weight+1] = max(items[item-1][1] + tab[item-1, weight+1 - items[item-1][2]], tab[item-1,weight+1])
			else # case 2: if it doesn't fit then just add the best previous
				tab[item,weight+1] = tab[item-1,weight+1]
			end
		end
	end
	return tab[numItems+1, Capacity+1]
end
#Vector{Tuple{Int64,Int64}}

Capacity = 5000
numItems = 5000
println(Capacity)
println(numItems)
items = Vector{Tuple{Int64,Int64}}(undef,numItems)
for i = 1:numItems
	items[i] = (rand(1:numItems), rand(1:Capacity))
	println(items[i][1], " ", items[i][2])
end
println(0)
ans = knpsk(Capacity, numItems, items)
println(ans)
