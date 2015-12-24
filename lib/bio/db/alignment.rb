# Extending Bio::DB::Alignment with mutation calling method
Bio::DB::Alignment.class_eval do
	# Aliases included for more intuitive naming when using for genomic alignments:
	alias_method :chr, :rname
	alias_method :opt, :tags # vice versa as opt is the "proper" sam name for the tag fields

	# def initialize(line)
	# 	@fields = line.split(/\s+/)
	# 	@qname, @flag, @rname, @pos, @mapq, cigar, @mrnm, @mpos, @isize, @seq, @qual = @fields
  #   @cigar = Bio::Alignment::CIGAR.new(cigar,@seq,source="sam")
	# 	@tags = @fields[11..-1]
	# 	@h = Hash.new{|k,v| @h[k] = ""}
	# 	@h[:unsplit_tags] = @tags.join(" ")
	# 	# Reconstruct the tag and assign a value. Assumes all tags are of the form X:Y:stuff
	# 	@tags.map!{|tag| tag.split(/:/,3)}
	# 	@tags.each{|a,b,c| @h[a + ":" + b] = c}
	# 	@tags = @h
	# 	@pos = @pos.to_i
	# 	@mapq = @mapq.to_i
	# 	@mrnm == "*" ? @mrnm = nil : @mrnm = @mrnm
	# 	@mpos == "0" ? @mpos = nil : @mpos = @mpos.to_i
	# 	@isize == "0" ? @isize = nil : @isize = @isize.to_i
	# 	@seq = Bio::Sequence::NA.new(@seq)
	# end

  # Extend getter methods for these as they will be set to objects by the mutations method
	# def cigar
	# 	@cigar.is_a? Bio::Alignment::CIGAR ? @cigar.string : @cigar
	# end
	#
	# def seq
	# 	@seq.is_a? Bio::Sequence::NA ? @seq.seq : @seq
	# end

  # Call mutations
  # Want to be able to give a length and offset - use this to generate appropriate sub CIGARs, subMDs & call
	def mutations offset=1, length=@seq.length, reference_pos=@pos-1, translation_start=1
		seq = Bio::Sequence::NA.new(@seq)
		cigar = Bio::Alignment::CIGAR.new(@cigar,seq,source="sam")
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
				  total += pair[1]
				when "I"
					mut = Bio::Mutation.new
					mut.type = :insertion
					mut.reference = nil
					mut.position = reference_pos
					mut.mutant = (insertions.length == 0) ? "N" : insertions.shift.upcase
					mut.seqname = @rname.to_s
          mutations << mut
			end
		end

    # Now substitutions & deletions - these need the MD tag
    sub_pos = mdz.report(/[sd]/)
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
    mutations.length > 0 ? mutations.sort{|x,y| x.position.to_i <=> y.position.to_i} : nil
  end
end