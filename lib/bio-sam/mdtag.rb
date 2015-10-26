#!/usr/bin/env ruby
# In the following examples, the guide is: CGTTCTTATGGAAGCCGGGA followed by AGG PAM.
#
# In all these cases, take the genome position and add the matched length from CIGAR string to get the genome position of the mutation.
#
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
#
# Methods for MD parsing:
# 	-cumulative reference length (eq. to matched_length in CIGAR).
# 	-organise into deletions & mismatches [no insertions in MD:Z]
#
require 'test/unit'
require 'rubygems'
require 'bio'
require '/Users/spettitt/work/scripts/cigar/iterate_pairs'

class Bio::Alignment::SAM
	# Field names from SAM specification
	attr_accessor :qname, :flag, :rname, :pos, :mapq, :cigar, :mrnm, :mpos, :isize, :seq, :qual, :tags
	# Aliases included for more intuitive naming when using for genomic alignments:
	alias_method :chr, :rname
	alias_method :opt, :tags # vice versa as opt is the "proper" sam name for the tag fields

	def initialize(line)
		@fields = line.split(/\s+/)
		@qname, @flag, @rname, @pos, @mapq, @cigar, @mrnm, @mpos, @isize, @seq, @qual = @fields
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

end

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
