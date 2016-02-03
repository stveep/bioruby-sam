class MutantAllele
  attr_accessor :mutations, :count, :example
  self << class
    attr_accessor: previous_lookups
  end
  @previous_lookups = {}

  def initialize (mutations: Bio::MutationArray.new, count: 0, example: nil)
  end

  # Returns JSON from Ensembl VEP
  def lookup species="human", ref_type=nil
    key = mutations.to_hgvs(ref_type)
    if key && MutantAllele.previous_lookups.keys.include? key
      MutantAllele.previous_lookups[key]
    else
      mutations.vep(species,ref_type)
    end
  end
end
