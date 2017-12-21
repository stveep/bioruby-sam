require 'helper'
require 'net/ping'
class MutationTest < Test::Unit::TestCase
  def test_mutation_initialisation
    params = {position: 1234,
              type: :deletion,
              reference: "ATGG",
              mutant: nil,
              seqname: 5}
    mut = Bio::Mutation.new(params)
    assert_equal 1234, mut.position
    assert_equal :deletion, mut.type
    assert_equal "ATGG", mut.reference
    assert_nil mut.mutant
    assert_equal 5, mut.seqname
  end

  def test_single_deletion_hgvs
    params = {position: 1234,
              type: :deletion,
              reference: "G",
              mutant: nil,
              seqname: 5}
    mut = Bio::Mutation.new(params)
    assert_equal "1234delG", mut.to_hgvs
    assert_equal "5:g.1234delG", mut.to_hgvs("g")
  end

  def test_deletion_hgvs
    params = {position: 1234,
              type: :deletion,
              reference: "ATGG",
              mutant: nil,
              seqname: 5}
    mut = Bio::Mutation.new(params)
    assert_equal "1234_1237delATGG", mut.to_hgvs
    assert_equal "5:g.1234_1237delATGG", mut.to_hgvs("g")
  end

  def test_substitution_hgvs
    params = {position: 1234,
              type: :substitution,
              reference: "A",
              mutant: "T",
              seqname: 5}
    mut = Bio::Mutation.new(params)
    assert_equal "1234A>T", mut.to_hgvs
    assert_equal "5:g.1234A>T", mut.to_hgvs("g")
  end

  def test_longer_substitution_hgvs
    params = {position: 1234,
              type: :substitution,
              reference: "ATG",
              mutant: "CCT",
              seqname: 5}
    mut = Bio::Mutation.new(params)
    assert_equal "1234_1236ATG>CCT", mut.to_hgvs
    assert_equal "5:c.1234_1236ATG>CCT", mut.to_hgvs("c")
  end

  def test_insertion_hgvs
    params = {position: 1234,
              type: :insertion,
              reference: nil,
              mutant: "T",
              seqname: 5}
    mut = Bio::Mutation.new(params)
    assert_equal "1234_1235insT", mut.to_hgvs
    assert_equal "5:g.1234_1235insT", mut.to_hgvs("g")
  end

  def test_vep
    skip_unless_ensembl_available do
      params = {position: 112839914,
                type: :deletion,
                reference: "ACCACC",
                mutant: nil,
                seqname: 5}
      mut = Bio::Mutation.new(params)
      mut_array = Bio::MutationArray.new([mut])
      result = mut_array.vep("human","g")
      assert_equal "inframe_deletion", result.first["most_severe_consequence"], "NB: Test written for GRCh38"
      assert_equal "ACCACC/-", result.first["allele_string"], "NB: Test written for GRCh38"
    end
  end
end
