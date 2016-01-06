# Parse a CIGAR string
# An example from Exonerate output. Ideally will also allow SAM file input to be used.
#   1 : CGGCTATGGGGTCGTGGGTCCCGCGTTG-CTCTGGGGCTCGGCACCCTGGGGCGGCACGGCCGT :  63
#       | | || | ||||||||||||||||||| |||||||||||||||||||||||||||||||||||
#   1 : CAG-TA-GTGGTCGTGGGTCCCGCGTTGTCTCTGGGGCTCGGCACCCTGGGGCGGCACGGCCGT :  62
#
# ref: CAGTAGTGGTCGTGGGTCCCGCGTTGTCTCTGG...
# cigar: SP-A12_D02_2015-01-16.seq 0 611 + SP-A3_ref 0 621 + 2514  M 3 I 1 M 2 I 1 M 21 D 1 M 306 D 9 M 89 I 1 M 126 D 1 M 24 D 1 M 8 I 1 M 6 D 1 M 5 D 1 M 17
# I are not counted in reference
# Regexp (from SAM specification) - but in exonerate the number comes first: ([0-9]+[MIDNSHP])+|\*

class Bio::Alignment::CIGAR
	include Bio::Alignment::IteratePairs
	@@regexps = {"exonerate" => /([MIDNSHP]{1})(\d+)/, "sam" => /(\d+)([MIDNSHP]{1})/}
	# Type of elements that count towards the reference length:
	# TODO: add full support for other elements S, H etc.
	@@reference = /[MD]/
	@@subexp = /([atgcAGCT]+)>([atgcAGTC]+)/
	attr_accessor :pairs, :reference

	def initialize(string,ref=nil,source="")
		# strip out whitespace
		string.gsub!(/\s+/,"")

		# Auto-detect source if not supplied
		if !(@@regexps.keys.include? source)
		  @@regexps.each do |k,v|
				# Look for match at start of string
				if m = string.match(v)
				  source = k if m.offset(0)[0] == 0
			  end
			end
			if source == ""
				raise "Source (e.g. 'exonerate', 'sam') not given and failed to auto-detect."
			end
    end
		# Make an array of pairs of of cigar elements:
		@pairs = string.scan(@@regexps[source])
		if source == "exonerate"
			@pairs.map!{|pair| [pair[0].to_s, pair[1].to_i]}
		else
			# Provision to have number and identifier the other way round
			@pairs.map!{|pair| [pair[1].to_s, pair[0].to_i]}
		end

		# Include reference sequence if provided
		@reference = ref
		# Check length of reference = sum(M+D)?
		#warn "Reference length is not equal to that implied by CIGAR string: #{@reference.length}, #{self.reference_length}." unless @reference.length == self.reference_length

	end

	# Given an offset in reference sequence and length, return an object corresponding to that subregion of the alignment
	def subalignment(offset,length,regexp=@@reference)
		new_array = iterate_pairs(@pairs,offset,length,regexp)
		# Return a CIGAR instance with just the new alignment
		new_string = new_array.join(" ")
		# -1 from offset as ruby string starts at zero
		new_cigar = Bio::Alignment::CIGAR.new(new_string,@reference[offset-1,length])
		new_cigar.remove_empty!
	end
	alias_method :slice, :subalignment

	# Given a CIGAR-based [not reference - use subalignment] offset and length, return a subregion
	def subcigar(offset,length)
		# No regexp - includes everything
		self.subalignment(offset,length,//)
	end

	def unmasked(offset,length)
		self.subalignment(offset,length,/[MDI]/)
	end

	def remove_small!(threshold=1)
		# Deletions convert to matches, insertions just remove
		deletions_to_matches(threshold)
		remove_small_nonmatches(threshold)
		self
	end

	def remove_empty!
		self.pairs.keep_if{|pair| pair[1] != 0 }
		self
	end

	def matched_length
		count_type("M")
	end

	def deleted_length
		count_type("D")
	end

	def inserted_length
		count_type("I")
	end

	def masked_length
		count_type(/[SH]/)
	end

	def reference_length
		count_type(@@reference)
	end

	def query_length
		count_type(/[MHSI]/)
	end

	# Output a representation of the query: replace deleted portions with "-", flag insertions with "*" or sim. Optionally provide the sequence (or symbols to use) of insertions, in order of appearence.
	# Should be able to accept an array
	# TODO: Add support for substitution highlighting (e.g lowercasing)
	def query(insertions=nil)
		if (insertions && (insertions.is_a? String))
			insertions = [insertions]
		end
		sequence = []
		total = 0
		@pairs.each do |pair|
			if pair[0].match("M")
				sequence << @reference[total..total+pair[1]-1].upcase
				total += pair[1]
			end
			if pair[0].match("I")
				if (insertions)
					insertion = insertions.shift.to_s
				else
					insertion = '['+pair[1].to_s+']'
				end
				sequence << insertion
			end
			if pair[0].match("D")
				pair[1].times{ sequence << "-" }
				total += pair[1]
			end
		end
		sequence.join("")
	end

	# Output hgnc variant format given reference position. Only deletions can be accurately annotated from the cigar string; insertions or wild type seqeunces return nil
	# NB mutation calling and annotation now implemented as extension to Bio::DB::Alignment (SAM)
	def hgnc(reference_pos=0,insertions=[],type="g",*subs)
		if insertions
			if insertions.is_a? String
				insertions = [insertions]
			end
		end
		first_match = true
		total = 0
		hgnc_format = []
		@pairs.each do |pair|
			case pair[0]
				when "M"
					#break if first_match == false
					reference_pos += pair[1]
					total += pair[1]
					first_match = false
				when "D"
					deleted_bases = @reference[total,pair[1]].upcase
					if (pair[1] == 1)
						string = (reference_pos + 1).to_s
					else
						string = (reference_pos + 1).to_s + "_" + (reference_pos + pair[1]).to_s
					end
					string = string + "del" + deleted_bases
					hgnc_format << string
					total += pair[1]
				when "I"
					inserted_bases = (insertions.length == 0) ? "N" : insertions.shift

					hgnc_format << (reference_pos).to_s + "_" + (reference_pos + 1).to_s + "ins" + inserted_bases.upcase
			end
		end
		# Use for substitutions, but could also pass any other annotation to include in here, as an array of strings
		subs = subs.first # >1 arguments discarded
		if subs
			if (subs.length > 0 && (subs.is_a? Array))
				hgnc_format = hgnc_format + subs
			end
		end
		if hgnc_format.length == 0
			nil
		elsif hgnc_format.length == 1
			type.to_s + "." + hgnc_format[0]
		else
			type.to_s + "." + "[" + hgnc_format.join(";") + "]"
		end

	end

  # TODO combine adjacent operations of the same type into a single pair
	def combine_adjacent

	end

	# Returns a hash (keyed by operation type) of three element arrays: the start positions on the reference of operations of the given type(s) and the length of the operation,
	# followed by query position (for e.g. retrieving inserted bases from SAM). A regexp can be used to specify multiple types e.g. /[ID]/.
	def positions(type)
		total = 0
		qtotal = 0
		hash = Hash.new{|h,k| h[k] = []}
		@pairs.each do |pair|
			if pair[0].match(type)
				hash[$&] << [total, pair[1], qtotal]
			end
			total += pair[1] if pair[0].match(@@reference)
			qtotal += pair[1]
		end
		hash
	end


	private

	def deletions_to_matches(threshold)
		self.pairs = @pairs.each{|pair| pair[0].sub!("D","M") if pair[1] <= threshold}
	end

	def remove_small_nonmatches(threshold)
		self.pairs = @pairs.keep_if{|pair| pair[0] == "M" || pair[1] > threshold}
	end

	def count_type(type)
		sum = 0
		@pairs.each do |pair|
			if pair[0].match(type)
				sum += pair[1]
			end
		end
		sum
	end

end #class
