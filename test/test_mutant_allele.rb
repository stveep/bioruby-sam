require 'helper'
class MutantAlleleTest < Test::Unit::TestCase
  def test_initialization
    ma = MutantAllele.new
    assert_equal 0, ma.count
  end

  def test_cached_lookup
    ma = MutantAllele.new
    m = Bio::Mutation.new
    m.position = 358
    m.reference = "TCC"
    m.seqname = "ENST00000366794"
    m.type = :deletion
    m.mutant = nil
    ma.mutations = Bio::MutationArray.new([m])
    MutantAllele.previous_lookups[ma.mutations.to_hgvs] = "mock previous lookup"
    assert_equal "mock previous lookup", ma.lookup
  end

end
