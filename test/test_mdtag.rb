require 'helper'
class MDZTest < Test::Unit::TestCase
	def test_match
		mdz = Bio::Alignment::SAM::MDZ.new("MD:Z:60^G13")
		assert_equal "60^G13", mdz.tag
	end

	def test_ref_length
		mdz_del = Bio::Alignment::SAM::MDZ.new("MD:Z:60^G13")
		assert_equal(mdz_del.ref_length, 74)
		mdz_sub = Bio::Alignment::SAM::MDZ.new("MD:Z:60G0A13")
		assert_equal(mdz_sub.ref_length, 75)
		mdz_none = Bio::Alignment::SAM::MDZ.new("MD:Z:150")
		assert_equal(mdz_none.ref_length, 150)

	end

	def test_pairs
		mdz = Bio::Alignment::SAM::MDZ.new("MD:Z:60^G13T")
		assert_equal [["m",60],["d","G"],["m",13],["s","T"]], mdz.pairs
		assert_equal [["m",60,0,0],["d","G",60,60],["m",13,61,60],["s","T",74,73]], mdz.cumulative
	end

	def test_reconstruct_tag
		mdz = Bio::Alignment::SAM::MDZ.new("MD:Z:60^G13T")
		assert_equal "60^G13T",  mdz.reconstruct_tag
	end

	def test_slice
		mdz = Bio::Alignment::SAM::MDZ.new("MD:Z:60^G13T")

		new_mdz = Bio::Alignment::SAM::MDZ.new("MD:Z:6^G3")
		assert_equal new_mdz.tag, mdz.slice(55,10).tag
	end

	def test_report
		mdz_del = Bio::Alignment::SAM::MDZ.new("MD:Z:60^G13")
		assert_equal([["d","G",60,60]], mdz_del.deletions)
		mdz_sub = Bio::Alignment::SAM::MDZ.new("MD:Z:60G0A13")
		assert_equal([["s","G",60,60],["s","A",61,61]], mdz_sub.substitutions)
		mdz_none = Bio::Alignment::SAM::MDZ.new("MD:Z:150")
		assert_equal([], mdz_none.report)
		mdz_del_sub = Bio::Alignment::SAM::MDZ.new("MD:Z:60^G5A13")
		assert_equal [["s","A",66,65]], mdz_del_sub.substitutions
	end
end
