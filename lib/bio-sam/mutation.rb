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
  def to_hgvs

  end

  # TODO Look up in Ensembl VEP, return JSON - see Ensembl REST API docs for details
  def vep

  end
end
