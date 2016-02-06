require 'helper'
class VepHgvsTest < Test::Unit::TestCase
  def test_vep_parsing
   json_obj = JSON.parse "[{\"assembly_name\":\"GRCh38\",\"end\":226392243,\"seq_region_name\":\"1\",\"transcript_consequences\":[{\"gene_id\":\"ENSG00000143799\",\"distance\":15,\"biotype\":\"protein_coding\",\"gene_symbol_source\":\"HGNC\",\"consequence_terms\":[\"downstream_gene_variant\"],\"strand\":-1,\"hgnc_id\":\"HGNC:270\",\"gene_symbol\":\"PARP1\",\"transcript_id\":\"ENST00000366792\",\"impact\":\"MODIFIER\"},{\"gene_id\":\"ENSG00000143799\",\"distance\":697,\"biotype\":\"protein_coding\",\"gene_symbol_source\":\"HGNC\",\"consequence_terms\":[\"downstream_gene_variant\"],\"strand\":-1,\"hgnc_id\":\"HGNC:270\",\"gene_symbol\":\"PARP1\",\"transcript_id\":\"ENST00000629232\",\"impact\":\"MODIFIER\"},{\"gene_id\":\"ENSG00000143799\",\"cdna_end\":227,\"biotype\":\"processed_transcript\",\"gene_symbol_source\":\"HGNC\",\"consequence_terms\":[\"non_coding_transcript_exon_variant\",\"non_coding_transcript_variant\"],\"strand\":-1,\"hgnc_id\":\"HGNC:270\",\"gene_symbol\":\"PARP1\",\"cdna_start\":225,\"transcript_id\":\"ENST00000469663\",\"impact\":\"MODIFIER\"},{\"cdna_end\":504,\"codons\":\"TCC/-\",\"protein_end\":120,\"strand\":-1,\"hgnc_id\":\"HGNC:270\",\"amino_acids\":\"S/-\",\"gene_symbol\":\"PARP1\",\"cdna_start\":502,\"transcript_id\":\"ENST00000366794\",\"cds_start\":358,\"gene_id\":\"ENSG00000143799\",\"protein_start\":120,\"biotype\":\"protein_coding\",\"gene_symbol_source\":\"HGNC\",\"cds_end\":360,\"consequence_terms\":[\"inframe_deletion\"],\"impact\":\"MODERATE\"}],\"strand\":-1,\"id\":\"ENST00000366794:c.358_360delTCC\",\"allele_string\":\"TCC/-\",\"most_severe_consequence\":\"inframe_deletion\",\"start\":226392241}]"
   #test result from hg18
   puts json_obj
   assert_equal [{"Allele"=>"TCC/-", "CDS position"=>358, "Protein start"=>120, "Mutation"=>"S/-", "Consequence"=>["inframe_deletion"]}], VepHgvs.consequences_for_transcript(json_obj,"ENST00000366794")

  end
end
