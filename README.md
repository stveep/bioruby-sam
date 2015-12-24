# bio-sam

[![Build Status](https://secure.travis-ci.org/stveep/bioruby-sam.png)](http://travis-ci.org/stveep/bioruby-sam)

Full description goes here

Note: this software is under active development!

## Installation

```sh
gem install bio-sam
```

## Usage

```ruby
require 'bio-sam-mutation'

  # NB must be tab-delimited
  insertion_and_deletion = Bio::DB::Alignment.new("I2M5K:00253:00406\t0\t5\t112839854\t70\t63M2I138M1D27M7S\t*\t0\t0\tCAGTGATCTTCCAGATAGCCCTGGACAAACCATGCCACCAAGCAGAAGTAAAACACCTCCACCATACCTCCTCAAACAGCTCAAACCAAGCGAGAAGTACCTAAAAATAAAGCACCTACTGCTGAAAAGAGAGAGAGTGGACCTAAGCAAGCTGCAGTAAATGCTGCAGTTCAGAGGGTCCAGGTTCTTCCAGATGCTGATACTTATTACATTTTGCCACGGAAAGTACTGCTGAGG\t@CDDDCCCCACACCCCCCCC?CCACCCC>A6;;;;7;;6;6;BC;;6;;;;;.;;>ADDA??;;;;;?CCACCCD>C??@CCCC>C@C;>?CCCC@C=::@:::::+:::/:CCC?>>>>CCCCDDD9CCCC@AB????=AB>??;?BB>@@@AA???CC<@@?????BB>??;;;B<BC;??8;6:A=@=@BBB;;;?<77//*08*088888*8=9=?B7;;4;??????????<\tPG:Z:novoalign\tAS:i:183\tUQ:i:183\tNM:i:3\tMD:Z:201^T27")

  insertion_and_deletion.mutations
  #=> [#<Bio::Mutation:0x007fa20b5b4fc8 @position=112839916, @type=:insertion, @reference=nil, @mutant="AT", @seqname="5">, #<Bio::Mutation:0x007fa20b5b4960 @position=112840055, @type=:deletion, @reference="T", @mutant=nil, @seqname="5">]

  insertion_and_deletion.mutations.first.to_hgvs("g")
  #=> "5:g.112839916_112839917insAT"

  puts YAML.dump(insertion_and_deletion.mutations.first.vep("human","g").first["transcript_consequences"].keep_if{|c| c["transcript_id"] == "ENST00000257430"})
  #---
  # - variant_allele: AT
  #   cdna_end: 4379
  #   codons: cca/ccATa
  #   protein_end: 1441
  #   strand: 1
  #   hgnc_id: HGNC:583
  #   amino_acids: P/PX
  #   gene_symbol: APC
  #   cdna_start: 4378
  #   transcript_id: ENST00000257430
  #   cds_start: 4322
  #   gene_id: ENSG00000134982
  #   protein_start: 1441
  #   biotype: protein_coding
  #   gene_symbol_source: HGNC
  #   cds_end: 4323
  #   consequence_terms:
  #   - frameshift_variant
  #   impact: HIGH
  # => nil

# E.g. of full request return
# http://rest.ensembl.org/documentation/info/vep_hgvs_get
insertion_and_deletion.mutations(112839854).first.vep("human","g")
# => [{"assembly_name"=>"GRCh38", "end"=>112839917, "seq_region_name"=>"5", "transcript_consequences"=>[{"gene_id"=>"ENSG00000134982", "distance"=>46, "variant_allele"=>"AT", "biotype"=>"nonsense_mediated_decay", "gene_symbol_source"=>"HGNC", "consequence_terms"=>["downstream_gene_variant"], "strand"=>1, "hgnc_id"=>"HGNC:583", "gene_symbol"=>"APC", "transcript_id"=>"ENST00000502371", "impact"=>"MODIFIER"}, {"variant_allele"=>"AT", "cdna_end"=>4380, "codons"=>"-/AT", "protein_end"=>1442, "strand"=>1, "hgnc_id"=>"HGNC:583", "amino_acids"=>"-/X", "gene_symbol"=>"APC", "cdna_start"=>4379, "transcript_id"=>"ENST00000257430", "cds_start"=>4323, "gene_id"=>"ENSG00000134982", "protein_start"=>1441, "biotype"=>"protein_coding", "gene_symbol_source"=>"HGNC", "cds_end"=>4324, "consequence_terms"=>["frameshift_variant"], "impact"=>"HIGH"}, {"gene_id"=>"ENSG00000134982", "distance"=>863, "variant_allele"=>"AT", "biotype"=>"protein_coding", "gene_symbol_source"=>"HGNC", "consequence_terms"=>["downstream_gene_variant"], "strand"=>1, "hgnc_id"=>"HGNC:583", "gene_symbol"=>"APC", "transcript_id"=>"ENST00000507379", "impact"=>"MODIFIER"}, {"variant_allele"=>"AT", "cdna_end"=>4481, "codons"=>"-/AT", "protein_end"=>1442, "strand"=>1, "hgnc_id"=>"HGNC:583", "amino_acids"=>"-/X", "gene_symbol"=>"APC", "cdna_start"=>4480, "transcript_id"=>"ENST00000508376", "cds_start"=>4323, "gene_id"=>"ENSG00000134982", "protein_start"=>1441, "biotype"=>"protein_coding", "gene_symbol_source"=>"HGNC", "cds_end"=>4324, "consequence_terms"=>["frameshift_variant"], "impact"=>"HIGH"}, {"gene_id"=>"ENSG00000134982", "distance"=>409, "variant_allele"=>"AT", "biotype"=>"protein_coding", "gene_symbol_source"=>"HGNC", "consequence_terms"=>["downstream_gene_variant"], "strand"=>1, "hgnc_id"=>"HGNC:583", "gene_symbol"=>"APC", "transcript_id"=>"ENST00000512211", "impact"=>"MODIFIER"}, {"gene_id"=>"ENSG00000134982", "variant_allele"=>"AT", "cdna_end"=>4569, "biotype"=>"nonsense_mediated_decay", "gene_symbol_source"=>"HGNC", "consequence_terms"=>["3_prime_UTR_variant", "NMD_transcript_variant"], "strand"=>1, "hgnc_id"=>"HGNC:583", "gene_symbol"=>"APC", "cdna_start"=>4568, "transcript_id"=>"ENST00000508624", "impact"=>"MODIFIER"}, {"gene_id"=>"ENSG00000258864", "variant_allele"=>"AT", "biotype"=>"nonsense_mediated_decay", "gene_symbol_source"=>"Clone_based_vega_gene", "consequence_terms"=>["intron_variant", "NMD_transcript_variant"], "strand"=>1, "gene_symbol"=>"CTC-554D6.1", "transcript_id"=>"ENST00000520401", "impact"=>"MODIFIER"}, {"gene_id"=>"ENSG00000134982", "distance"=>2195, "variant_allele"=>"AT", "biotype"=>"protein_coding", "gene_symbol_source"=>"HGNC", "consequence_terms"=>["downstream_gene_variant"], "strand"=>1, "hgnc_id"=>"HGNC:583", "gene_symbol"=>"APC", "transcript_id"=>"ENST00000504915", "impact"=>"MODIFIER"}], "strand"=>1, "id"=>"5:g.112839917_112839918insAT", "allele_string"=>"-/AT", "most_severe_consequence"=>"frameshift_variant", "start"=>112839918}]



```

The API doc is online. For more code examples see the test files in
the source tree.

## Project home page

Information on the source tree, documentation, examples, issues and
how to contribute, see

  http://github.com/stveep/bioruby-sam

The BioRuby community is on IRC server: irc.freenode.org, channel: #bioruby.

## Cite

If you use this software, please cite one of

* [BioRuby: bioinformatics software for the Ruby programming language](http://dx.doi.org/10.1093/bioinformatics/btq475)
* [Biogem: an effective tool-based approach for scaling up open source software development in bioinformatics](http://dx.doi.org/10.1093/bioinformatics/bts080)

## Biogems.info

This Biogem is published at (http://biogems.info/index.html#bio-sam)

## Copyright

Copyright (c) 2015 stveep. See LICENSE.txt for further details.
