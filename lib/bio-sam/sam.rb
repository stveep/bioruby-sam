# Insertions only appear in the cigar string; but the inserted base can only be found from the read  e.g.
#DKNQZ:00014:00145       0       5       112767204       37      58M1I39M        *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCCGGGAAGGATCTGTATCAAGCCGTTCTGGAGAGTGCAG      CBDC??<;D1;?=?C@CC;;;;/;;;;;@C;;;;;-;;444*4@79@.22499:<@<:::/=D<D@D9=;;;;;;;44.2274;)-.-)--332226=      XT:A:U  NM:i:1  X0:i:1  X1:i:0  XM:i:0  XO:i:1  XG:i:1  MD:Z:97
# in the above case the inserted bases would be sequence_string[total matched length .. total matched length+insertion length] i.e. sequence_string[58..59] => "C"
#
# Deletions will be in both - they are ^N in the MD:Z tag, e.g
# DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U  NM:i:3  X0:i:1  X1:i:0  XM:i:3  XO:i:1  XG:i:1  MD:Z:60^G13
# To get the base, need to parse the MD:Z tag and take the [ATGC] characters after the ^.
#
# In the case of a mismatch, the CIGAR string will be M (As this is made for alignment, not variant calling), whereas the reference base will appear in the MD:Z tag - e.g. here when CCgGG in the reference becomes CCCGG:
# DKNQZ:00167:00152       0       5       112767204       25      157M    *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCCGGAAGGATCTGTATCAAGCCGTTCTGGAGAGTGCAGTCCTGTTCCTATGGGTTCATTTCCAAGAAGAGGGTTTGTAAATGGAAGCAGAGGCTGAGG   ??CCCCACD>CC>CC>BB;///(00/--3--2222*222::=:C65882>>A::5:5:AC;;6:66).-33322529.29/22/9:2A522CC;:?CCEACCD@@@@DCDG@FA??<<;?@DDFDCC@@@=@D<@@CC>CC>?>>AC@;>6;A@;;7   XT:A:U  NM:i:7  X0:i:1  X1:i:0  XM:i:7  XO:i:0  XG:i:0  MD:Z:60G89A0A0A1T0A0C0
#
# ===>  Need the MD, CIGAR and reference sequence
class Bio::Alignment::SAM
	# Field names from SAM specification
	attr_accessor :qname, :flag, :rname, :pos, :mapq, :cigar, :mrnm, :mpos, :isize, :seq, :qual, :tags
	# Aliases included for more intuitive naming when using for genomic alignments:
	alias_method :chr, :rname
	alias_method :opt, :tags # vice versa as opt is the "proper" sam name for the tag fields

	def initialize(line)
		@fields = line.split(/\s+/)
		@qname, @flag, @rname, @pos, @mapq, cigar, @mrnm, @mpos, @isize, @seq, @qual = @fields
    @cigar = Bio::Alignment::CIGAR.new(cigar,@seq,source="sam")
		@tags = @fields[11..-1]
		@h = Hash.new{|k,v| @h[k] = ""}
		@h[:unsplit_tags] = @tags.join(" ")
		# Reconstruct the tag and assign a value. Assumes all tags are of the form X:Y:stuff
		@tags.map!{|tag| tag.split(/:/,3)}
		@tags.each{|a,b,c| @h[a + ":" + b] = c}
		@tags = @h
		@pos = @pos.to_i
		@mapq = @mapq.to_i
		@mrnm == "*" ? @mrnm = nil : @mrnm = @mrnm
		@mpos == "0" ? @mpos = nil : @mpos = @mpos.to_i
		@isize == "0" ? @isize = nil : @isize = @isize.to_i
		@seq = Bio::Sequence::NA.new(@seq)
	end

  # Call mutations
  # Want to be able to give a length and offset - use this to generate appropriate sub CIGARs, subMDs & call
	def mutations reference_pos=0, offset=1, length=@seq.length, translation_start=1
    # Generate subalignments from the CIGAR and MD:Z
    subcigar = @cigar.subalignment(offset,length)
    mdz = Bio::Alignment::SAM::MDZ.new(@tags["MD:Z"])
		mdz = mdz.slice(offset,length)
    # Get inserted bases from the read sequence, only within the region of interest
    insertions = []
    insertion_positions = subcigar.positions(/I/)
    unless insertion_positions.empty?
      insertion_positions["I"].each do |ins|
        # Sam.seq returns a Sequence::NA object
        # Need a -1 as ruby counts characters
        # Use ins[2] to retrieve the base as this is the position on query. ins[1] is position on reference, used to annotate position.
        i = @seq.seq[(offset+ins[2]-1),ins[1]]
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
				preceding = @cigar.subalignment(0,offset-1)
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
            mut.mutant = @seq[read_position,p[1].length].upcase
            mutations << mut
          when "d"
            mut = Bio::Mutation.new
  					mut.type = :deletion
  					mut.reference = p[1].upcase
  					mut.position = substart+p[2] + 1
  					mut.mutant = nil
  					mutations << mut
        end
      end
    end
    mutations.length > 0 ? mutations.sort{|x,y| x.position.to_i <=> y.position.to_i} : nil
  end

end
