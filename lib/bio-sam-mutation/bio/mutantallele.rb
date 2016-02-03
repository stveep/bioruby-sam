class MutantAllele
  attr_accessor :mutations, :count, :example
  class << self
    attr_accessor :previous_lookups
  end
  @previous_lookups = {}

  def initialize (mutations: nil, count: 0, example: nil)
    @mutations = mutations
    @count = count
    @example = example
  end

  # Returns JSON from Ensembl VEP
  def lookup species="human", ref_type=nil
    key = mutations.to_hgvs(ref_type)
    if key && (MutantAllele.previous_lookups.keys.include? key)
      MutantAllele.previous_lookups[key]
    else
      mutations.vep(species,ref_type)
    end
  end
end
