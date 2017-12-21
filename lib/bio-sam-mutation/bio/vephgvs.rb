module VepHgvs
  require 'bio-ensembl-rest'
  def vep(species="human",reference_type=nil)
    return unless self.first.to_hgvs
    if reference_type.nil?
      reference_type = self.first.seqname.match(/^ENS/) ? "c" : "g"
    end
    EnsemblRest.connect_db
    JSON.parse(EnsemblRest::Variation.vep_hgvs(species,self.to_hgvs(reference_type)))
  end

  def self.consequences_for_transcript (json_object,tscript)
    if json_object.length > 0
      if json_object.first["transcript_consequences"]
        consequences = json_object.first["transcript_consequences"]
        consequences.keep_if{|a| a["transcript_id"] == tscript}
        consequences.map{|a| {"Allele" => json_object.first["allele_string"], "CDS position" => a["cds_start"], "Protein start" => a["protein_start"], "Mutation" => a["amino_acids"], "Consequence" => a["consequence_terms"]}}
      end
    end
  end
end
