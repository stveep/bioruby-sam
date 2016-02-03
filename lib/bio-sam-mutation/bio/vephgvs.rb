module VepHgvs
  require 'bio-ensembl-rest'
  def vep(species="human",reference_type=nil)
    if reference_type.nil?
      reference_type = @seqname.match(/^ENS/) ? "c" : "g"
    end
    EnsemblRest.connect_db
    JSON.parse(EnsemblRest::Variation.vep_hgvs(species,self.to_hgvs(reference_type)))
  end

  def self.consequences_for_transcript (json_object,tscript)
    if json_object.length > 0
      if json_object.first["transcript_consequences"]
        consequences = json_object.first["transcript_consequences"]
        cons_of_interest = consequences.keep_if{|a| a["transcript_id"] == tscript}
        cons_of_interest.map!{|a| { "CDS position" => a["cds_start"], "Protein start" => a["protein_start"], "Mutation" => a["amino_acids"], "Consequence" => a["consequence_terms"]}}
      end
    end
  end
end
