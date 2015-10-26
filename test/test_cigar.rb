require 'helper'
class CigarTest < Test::Unit::TestCase
	def test_partition
		s = Bio::Alignment::CIGAR.new("M 1 D 2 I 3","GAT",source="exonerate")
		assert_equal([["M",1],["D",2],["I",3]], s.pairs)
		assert_equal(3, s.pairs.length)
		assert_equal(1, s.matched_length)
	end

	def test_auto_detect
		exonerate  = Bio::Alignment::CIGAR.new("M 1 D 2 I 3","GAT")
		sam =  Bio::Alignment::CIGAR.new("2M5D4M","ATCCTCCGGAA")
		assert_raise {Bio::Alignment::CIGAR.new("rubbish","GATC")}
		assert_nothing_raised {Bio::Alignment::CIGAR.new("M 1 D 2 I 3","GAT")}
		assert_equal [["M",1],["D",2],["I",3]], exonerate.pairs
		assert_equal [["M",2],["D",5],["M",4]], sam.pairs
	end

	def test_subalignment
		m = Bio::Alignment::CIGAR.new("M 35","GGATCGATCGATCGATCGATCGATCGATCGATCGA")
		msa = m.subalignment(10,12)
		assert_equal "GATCGATCGATC", msa.reference
		assert_equal [["M", 12]], msa.pairs

		s = Bio::Alignment::CIGAR.new("M 16 I 1 M 325 D 3 M 4 D 2 M 178 I 1 M 17","TTTAATTGCATTTAATTGCATTTAATTGCATTAATTGCATGGTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTA")
		sa = s.subalignment(340,7)
		# M 16 + M 325 [I not included] = 341
		# Start at 340 = include 340 and 341, i.e. M 2
		# Total of all M+D in subalignment should be 7
		assert_equal([["M", 2],["D",3],["M",2]],sa.pairs)
		assert_equal("CATTTAA",sa.reference)
		# Case where subalignment starts at the end of a pair (M 325 in test case)
		sa2 = s.subalignment(341,7)
		assert_equal([["M", 1],["D",3],["M",3]],sa2.pairs)
	end

	def test_subcigar
		s = Bio::Alignment::CIGAR.new("M 16 I 1 M 325 D 3 M 4 D 2 M 178 I 1 M 17","TTTAATTGCATTTAATTGCATTTAATTGCATTAATTGCATGGTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTAATTGCATTTA")
		soft = Bio::Alignment::CIGAR.new("1S5M6D15M","CCCTGACGTTGAGGTGGATGGGTTCTCTGAGCTTCGG",source="sam")
		assert_equal [["M", 2],["I",1],["M",4]], s.subcigar(15,7).pairs
		assert_equal [["M",16],["I",1]], s.subcigar(1,17).pairs
		assert_equal [["S",1],["M",2]], soft.subcigar(1,3).pairs
	end

	def test_small
		s = Bio::Alignment::CIGAR.new("M 20 D 1 I 1 M 1 D 8 I 2 M 10","GGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT")
		assert_equal([["M",20],["M",1],["M",1],["D",8],["I",2],["M",10]],s.remove_small!.pairs)
		# Could improve by combining adjacent pairs of same type.
	end

	def test_small_deletions
		d = Bio::Alignment::CIGAR.new("M 20 D 1 M 10 D 4","GGATCGATCGATCGATCGATCGATCGATCGATCGA")
		i = Bio::Alignment::CIGAR.new("M 20 I 1 M 10 D 4","GGATCGATCGATCGATCGATCGATCGATCGATCG")
		assert_equal([["M",20],["M",1],["M",10],["D",4]],d.remove_small!.pairs)
		assert_equal([["M",20],["M",10],["D",4]],i.remove_small!.pairs)

	end

	def test_lengths
		m = Bio::Alignment::CIGAR.new("M 1 D 2 M 3 I 4 M 5","AAAATAAAATA")
		d = Bio::Alignment::CIGAR.new("M 1 D 2 M 3 I 4 D 5","AAAATAAAATA")
		i = Bio::Alignment::CIGAR.new("M 1 I 2 M 3 I 4 D 5","AAAATAAAA")
		assert_equal(9,m.matched_length)
		assert_equal(7,d.deleted_length)
		assert_equal(6,i.inserted_length)
		assert_equal(11,m.reference_length)
		assert_equal(11,d.reference_length)
		assert_equal(9,i.reference_length)
	end

  def test_masking
		s1 = Bio::Alignment::CIGAR.new("1S10M","AAAATAAAA",source="sam")
		s2 = Bio::Alignment::CIGAR.new("10M5D","AAAATAAAA",source="sam")
		h = Bio::Alignment::CIGAR.new("25H1I10M","AAAATAAAA",source="sam")
		assert_equal 1, s1.masked_length
		assert_equal 0, s2.masked_length
		assert_equal 25, h.masked_length
	end

	def test_query_length
		m = Bio::Alignment::CIGAR.new("M 1 D 2 M 3 I 4 M 5","AAAATAAAATA")
		d = Bio::Alignment::CIGAR.new("M 1 D 2 M 3 I 4 D 5","AAAATAAAATA")
		i = Bio::Alignment::CIGAR.new("M 1 I 2 M 3 I 4 D 5","AAAATAAAA")
		s1 = Bio::Alignment::CIGAR.new("1S10M","AAAATAAAA",source="sam")
		s = Bio::Alignment::CIGAR.new("1S5M2D2M1I5M","AAAATAAAA",source="sam")
		assert_equal 13, m.query_length
		assert_equal 8, d.query_length
		assert_equal 10, i.query_length
		assert_equal 11, s1.query_length
		assert_equal 14, s.query_length
	end

	def test_query_output
		m = Bio::Alignment::CIGAR.new("M 11","AAAATAAAATA")
		lc = Bio::Alignment::CIGAR.new("M 11","aaaataaaata")
    d = Bio::Alignment::CIGAR.new("M 5 D 2 M 4","AAAATAAAATA")
    i = Bio::Alignment::CIGAR.new("M 5 I 10 M 6","AAAATAAAATA")
    i2 = Bio::Alignment::CIGAR.new("M 5 I 2 M 3 I 4 M 6","AAAATGCCAAAATA")
		assert_equal("AAAATAAAATA",m.query)
		assert_equal("AAAATAAAATA",lc.query)
		assert_equal("AAAAT--AATA",d.query)
		assert_equal("AAAAT[10]AAAATA",i.query)
		assert_equal("AAAAT-insertion-AAAATA",i.query("-insertion-"))
		assert_equal("AAAATttGCCataaAAAATA",i2.query(["tt","ataa"]))
	end

	def test_positions
    i = Bio::Alignment::CIGAR.new("M 5 I 10 M 6","AAAATAAAATA")
		assert_equal({"I" => [[5, 10, 5]]}, i.positions(/I/))
    d = Bio::Alignment::CIGAR.new("M 5 D 2 M 4","AAAATAAAATA")
		assert_equal({"M" => [[0,5,0],[7,4,7]], "D" => [[5,2,5]]},d.positions(/[MD]/))
	end
	def test_positions_on_query
		i = Bio::Alignment::CIGAR.new("M 5 I 10 D 3 M 6","AAAATAAAATA")
		assert_equal({"D" => [[5, 3, 15]]}, i.positions(/D/))
	end

	def test_hgnc
		d = Bio::Alignment::CIGAR.new("M 5 D 2 M 4","AAAATAAAATA")
		i = Bio::Alignment::CIGAR.new("M 5 I 10 M 6","AAAATAAAATA")
		i2 = Bio::Alignment::CIGAR.new("M 5 I 2 M 3 D 1 M 6","AAAATGCCGAAAATA")
		dsam = Bio::Alignment::CIGAR.new("155M6D15M","CCCTGACGTTGAGGTGGATGGGTTCTCTGAGCTTCGGTGGGATGACCAGCAGAAAGTCAAGAAGACAGCGGAAGCTGGAGGAGTGACAGGCAAAGGCCAGGATGGAATTGGTAGCAAGGCAGAGAAGACTCTGGGTGACTTTGCAGCAGAGTATGCCAAGTCCAACAGAAGTACGT",source="sam")
		assert_equal "g.6_7delAA", d.hgnc
		assert_equal "g.25_26insGATCGATCGA", i.hgnc(20,"GATCGATCGA")
		assert_equal "g.[5_6insGA;9delG]", i2.hgnc(0,"GA")
		assert_equal "g.156_161delCCAAGT", dsam.hgnc(0)
	end

	def test_hgnc_from_subalignment
		align = Bio::Alignment::CIGAR.new("M 5 I 2 M 3 D 1 M 6","AAAATGCCGAAAATA")
		sa = align.subalignment(6,7)
	  assert_equal([["M", 3],["D",1],["M",3]],sa.pairs)
		assert_equal("GCCGAAA",sa.reference)
		assert_equal "g.4delG", sa.hgnc
	end

	def test_substitutions
		d = Bio::Alignment::CIGAR.new("M 5 D 2 M 4","AAAATAAAATA")
		s = Bio::Alignment::CIGAR.new("M 11","AAAATAAAATA")
		assert_equal "c.5T>A", s.hgnc(0,[],"c",["5T>A"])
		assert_equal "g.[6_7delAA;1A>T]", d.hgnc(0,[],"g",["1A>T"])

	end

end
