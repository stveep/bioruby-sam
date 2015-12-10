# Inherit from hash for easy json serialisation
class Bio::Mutation < Hash
  attr_accessor :position, :type, :reference, :mutant, :seqname
  def initialize params={position: 1,type: :uninitialized, reference: nil, mutant: nil, seqname:nil}
    @position = params[:position]
    @type = params[:type]
    @reference = params[:reference]
    @mutant = params[:mutant]
    @seqname = params[:seqname]
  end

  # TODO Annotate in HGVS format, return String
  # http://www.hgvs.org/mutnomen/recs.html
  # This gives just the annotation. To convert to a full allele description, needs to be combined
  # with e.g. g. for genomic:  g. - can supply this "g", "c" as type to annotate a single mutation directly
  # for compound mutants, need to join an array of annotations e.g. 1:g.[213456A>C;213460_213461delTG]
  def to_hgvs(reference_type=nil)
    if reference_type
      hgvs_arr = [@seqname,":",reference_type,".",@position.to_s]
    else
      hgvs_arr = [@position.to_s]
    end

    case @type
      when :deletion
        if @reference.length == 1
          hgvs_arr << "del"+@reference
        else
          hgvs_arr = hgvs_arr + ["_",
                                (@position.to_i+@reference.length-1).to_s,
                                "del",
                                @reference]
        end
        hgvs_arr.join

      when :substitution
        if @reference.length > 1
          hgvs_arr = hgvs_arr + ["_",
                                (@position.to_i+@reference.length-1).to_s]
        end
        hgvs_arr << @reference+">"+@mutant
        hgvs_arr.join

      when :insertion
        hgvs_arr << "_" + (@position.to_i+1).to_s
        hgvs_arr << "ins"+@mutant
        hgvs_arr.join
        # TODO - distinguish duplications from insertions? Needs further input from ref.
    end
  end

  # TODO Look up in Ensembl VEP, return JSON - see Ensembl REST API docs for details
  def vep

  end
end
