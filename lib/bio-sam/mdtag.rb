class Bio::Alignment::SAM::MDZ
	include Bio::Alignment::IteratePairs
	attr_accessor :tag, :pairs, :cumulative
	@@regexp = /MD:Z:([\w^]+)/
	@@format = /[\w^]+/
	@@splitter = /(?<match>\d+)|(?<substitution>[GATC]+)|\^(?<deletion>[GATC]+)/
	# Operations that consume reference seqeunce:
	@@reference = /[msd]/
	def initialize(string)
		if string.match(@@regexp)
			@tag = $~[1]
		elsif string.match(@@format)
			#Assume tag given without MD:Z: leader
			@tag = string
		else
			raise "Tag not of expected format."
		end
		# Splits the string into operations using the splitter regexp class variable, returns array of two-element arrays describing operations
		spl = @tag.scan(@@splitter)
		# Returns an array of matches [match,substition,deletion]
		# Although regexp captures are named, these don't get included automatically with scan as it doesn't return MatchData objects.
		spl.map! do |a|
			array = [["m", a[0]],["s", a[1]],["d", a[2]]]
			# Only one of these will be non-nil
			array.keep_if{|i| i[1]}
			array.map!{|i| if i[0] == "m" then i[1] = i[1].to_i end; i}
			array[0]
		end
		@pairs = spl

		@cumulative = []
		cumulative_length = 0
		read_length = 0
		@pairs.each do |q|
			p = q.dup
			case p[0]
				when "m"
					len = p[1]
					rlen = p[1]
				when "s"
					len = p[1].length
					rlen = p[1].length
				when "d"
					len = p[1].length
					# Deleted bases don't appear in the read, so don't count to the length
					rlen = 0
			end
			# third element in each array will be the total preceding length on the reference, i.e. the position of the operation.
			# fourth element is similar for the read.
			@cumulative << p.dup.push(cumulative_length).push(read_length)
			cumulative_length += len
			read_length += rlen
		end
	end

	def deletions
		report(/d/)
	end

	def substitutions
		report(/s/)
	end

	# Report the positions of given events
	def report(regexp=/[sd]/)
		to_return = []
		@cumulative.each do |p|
			if p[0] =~ regexp
				to_return << p
			end
		end
		to_return
	end

	# Reconstruct a MD:Z tag from the pairs array
	def reconstruct_tag(array=@pairs)
		new_tag = []
		array.each do |p|
			case p[0]
				when "m"
					string = p[1].to_s
				when "s"
					string = p[1]
				when "d"
					string = "^"+p[1]
			end
			new_tag << string
		end
		new_tag.join("")
	end



	# Sums the total length of the reference sequence represented by the MD:Z tag (or part of)
	def ref_length
		#Need the sum of all "movement" operations (i.e. numbers) as well as any substituted bases (count 1 each)
		if @tag =~ /^\d+$/
			@tag.to_i
		else
			temp_tag = @tag.dup
			temp_tag.gsub!(/\^/,"")  # Deletions need to be counted - sub the caret character out and count the remaining base characters
			movements = temp_tag.split(/[GATC]+/).map(&:to_i).reduce(:+) # Sum numbers
			deletions = temp_tag.split(/\d+/).map(&:length).reduce(:+) # Sum number of base chars
			movements + deletions
		end
	end
	# Given an offset in reference sequence and length, return an object corresponding to that subregion of the alignment
	def slice(offset,length)
		new_array = iterate_pairs(@pairs,offset,length,@@reference)
		# Return a MDZ instance with just the new alignment
		new_tag = reconstruct_tag(new_array)
		Bio::Alignment::SAM::MDZ.new(new_tag)
	end

end


# DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U  NM:i:3  X0:i:1  X1:i:0  XM:i:3  XO:i:1  XG:i:1  MD:Z:60^G13
