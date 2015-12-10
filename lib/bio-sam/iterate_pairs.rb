module Bio::Alignment::IteratePairs
	private
	#Mixin to iterate through ordered paired [operation, value] data and take subsets - e.g. broken down CIGAR and MD:Z tags. Can set a regexp for the operation; default matches everything
	def iterate_pairs(pairs,offset,length,regexp = //)
		offset = offset.to_i
		length - length.to_i
		total = 0
		new_array = []
		first = true
		pairs.each do |pair|
			new_pair = pair.dup
			if pair[1].is_a? String
				pairlength = pair[1].length
			elsif pair[1].is_a? Integer
				pairlength = pair[1]
			else
				raise "Value for operation must be a string or integer"
			end
			# Only count pairs where first element matches a regexp.
			# e.g. for CIGAR:
			# ref M + ref D = ref length
			# query M + query I = query length
			if pair[0].match(regexp)
				total += pairlength
			end
			# Just keep going until we get to the start of the subalignment
			if total < offset
				next
			end
			# If the offset is partway through a pair, need to split it up.
			if first
				# adjust the number in this pair; it will be added below
				if pair[1].is_a? String
					new_pair[1] = new_pair[1][total-offset..pairlength]
				else
					new_pair[1] = total - offset + 1  # Add one for bases
				end
			end

			# Once we are at/beyond the end of the desired region:
			if total >= offset + length
				# Special case where the whole subalignment is contained within one cigar element:
				if first
					if pair[1].is_a? String
						new_pair[1] = new_pair[1][total-offset-pairlength,length]
					else
						new_pair[1] = length
					end
					new_array << new_pair
					break
				end
			# Adding the last part of the alignment
				previous_total = total - pairlength
				if pair[1].is_a? Integer
					new_pair[1] = offset + length - previous_total - 1 #-1 extra for base arithmetic
				else
					new_pair[1] = new_pair[1][0..offset + length - previous_total]
				end
				new_array << new_pair
				break
			end
			first = false
			new_array << new_pair
		end
		new_array
	end
end
