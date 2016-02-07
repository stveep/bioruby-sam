# Extending Bio::DB::Alignment with mutation calling method
Bio::DB::Alignment.class_eval do
	# Aliases included for more intuitive naming when using for genomic alignments:
	alias_method :chr, :rname
	alias_method :opt, :tags # vice versa as opt is the "proper" sam name for the tag fields
	attr_accessor :cigar_obj

	def add_tag!(new_tag)
		if new_tag.is_a? String
			new_tag = Bio::DB::Tag.new(new_tag)
			# new_tag_obj.set(new_tag)
			# new_tag = new_tag_obj
		else
			raise "Tag not recognised - pass a string or Bio::DB::Tag object" unless new_tag.is_a? Bio::DB::Tag
		end
		@tags[new_tag.tag] = new_tag

		regenerate_string
	end

	def add_tag(tag)
		dup.add_tag!(tag)
	end

	# Output a representation of the query sequence
	def query offset=1, length=@seq.length, reference_pos=@pos-1, ins_chr="_"
	  mutations = self.mutations(offset,length)
		cigar = Bio::Alignment::CIGAR.new(@cigar,seq,source="sam")
		preceding = cigar.subalignment(0,offset-1)
		preceding_diff = preceding.query_length-(offset-1)
		pointer = preceding.query_length
		output = []
		deletions = 0
		insertions = 0
		if mutations
			mutations.each do |mut|
				mut.position = mut.position + insertions - deletions + preceding_diff
				case mut.type
					when :deletion
						# position for deletion is the first deleted base
						fillin = mut.position-1-reference_pos-1
				    output << @seq[pointer..fillin] if fillin > pointer
						mut.reference.length.times{ output << "-" }
						pointer += fillin - pointer + 1
						deletions += mut.reference.length
					when :insertion
						# position for insertion is the base we want
						fillin = mut.position-reference_pos-1
						output << @seq[pointer..fillin] if fillin > pointer
						output << ins_chr + mut.mutant.downcase + ins_chr
						pointer += fillin - pointer + 1 + mut.mutant.length
						insertions += mut.mutant.length
					when :substitution
						# position for substitution is the first subbed base
						fillin = mut.position-1-reference_pos-1
				    output << @seq[pointer..fillin] if fillin > pointer
						output << mut.mutant.downcase
						pointer += fillin - pointer + 1 + mut.mutant.length
				end
			end
		end
		# Remaining sequence
		if offset + length > pointer
      output << @seq[pointer..offset-1+length-1-deletions+insertions+preceding_diff]
		end
		output.join
	end

  # Call mutations
  # Want to be able to give a length and offset - use this to generate appropriate sub CIGARs, subMDs & call
	def mutations offset=1, length=nil, translation_start=1
		return nil if @query_unmapped
		cigar = Bio::Alignment::CIGAR.new(@cigar,seq,source="sam")
		length ||= cigar.reference_length - offset
		return nil if offset+length > cigar.reference_length
		seq = Bio::Sequence::NA.new(@seq)
		@cigar_obj = cigar
    # Generate subalignments from the CIGAR and MD:Z
    subcigar = cigar.subalignment(offset,length)
    mdz = Bio::DB::Tag::MD.new(@tags["MD"].value)
		mdz = mdz.slice(offset,length)
    # Get inserted bases from the read sequence, only within the region of interest
    insertions = []
    insertion_positions = subcigar.positions(/I/)
    unless insertion_positions.empty?
      insertion_positions["I"].each do |ins|
        # Sam.seq returns a Sequence::NA object
        # Need a -1 as ruby counts characters
        # Use ins[2] to retrieve the base as this is the position on query. ins[1] is position on reference, used to annotate position.
        i = seq.seq[(offset+ins[2]-1),ins[1]]
        insertions << i
      end
    end

    first_match = true
		total = 0
		mutations = []
		reference_pos = @pos - 1
		subcigar.pairs.each do |pair|
			case pair[0]
				when "M"
					#break if first_match == false
					reference_pos += pair[1]
					total += pair[1]
					first_match = false
        # Call deletions using the MD:Z tag - avoid need to supply reference seq.
				when "D"
				# Deletions are called below but still need to count here
				  reference_pos += pair[1]
				when "I"
					mut = Bio::Mutation.new
					mut.type = :insertion
					mut.reference = nil
					mut.position = reference_pos + offset - translation_start
					bases = insertions.shift
					mut.mutant = bases ? bases.upcase : "N"
					mut.seqname = @rname.to_s
          mutations << mut
			end
		end

    # Now substitutions & deletions - these need the MD tag
    sub_pos = mdz.report(/[sd]/)
		previous_sub_position = 0
		unless sub_pos.empty?
			sub_pos.each do |p|
				# Reference base is in the MD:Z tag (p[1] here), for the actual base need to go to the read
				# p[3] is the length of operations preceding the substitution on the read, p[2] on the reference.
				# p[2] and p[3] are defined on the subalignment, so should add them onto the preceding.
				# Need to add in any inserted bases from the CIGAR string using query_length
				preceding = cigar.subalignment(0,offset-1)
				# Masked length is not included in the MD:Z string so need to add it
				read_position = preceding.query_length+preceding.masked_length+p[3]
        # This is the adjustment needed to get the correct annotation:
        substart = @pos + offset - translation_start - 1
        case p[0]
          when "s"
            mut = Bio::Mutation.new
            mut.type = :substitution
            mut.position = substart+p[2] + 1
            mut.reference = p[1].upcase
            mut.mutant = seq[read_position,p[1].length].upcase
						mut.seqname = @rname.to_s
            mutations << mut
          when "d"
            mut = Bio::Mutation.new
  					mut.type = :deletion
  					mut.reference = p[1].upcase
  					mut.position = substart+p[2] + 1
  					mut.mutant = nil
						mut.seqname = @rname.to_s
  					mutations << mut
        end
      end
    end
    # mutations.length > 0 ? mutations.sort{|x,y| x.position.to_i <=> y.position.to_i} : nil
    mutations.length > 0 ? Bio::MutationArray.new(mutations.sort) : nil

  end

	def regenerate_string
		tags_string = @tags.map{|k,v| [v.tag, v.type, v.value].join(":") }
		self.sam_string = [@qname,
	        @flag,
	        @rname,
	        @pos,
	        @mapq,
	        @cigar,
	        @mrnm,
	        @mpos,
	        @isize,
	        @seq,
	        @qual,
	        tags_string].join("\t")
	end
end
