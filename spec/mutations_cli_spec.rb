describe MutationsCLI do
  describe "setting configuration" do
    it "sets up defaults for missing values in the config" do
      expect(MutationsCLI.set_defaults({file: "filename"}))
        .to eq({start: ".", 
               offset: 1, 
               length: 100, 
               translation_start: 1, 
               file:/filename/ })
    end

    it "doesn't replace values that have already been set, and doesn't set a filename regexp unless one was given" do
      expect(MutationsCLI.set_defaults({start: "TAG", length: 400}))
        .to eq({start: "TAG", 
               offset: 1, 
               length: 400, 
               translation_start: 1
              })
    end

    it "organises the configuration by product" do
      expect(MutationsCLI.construct_products({ products: {"prod1" => {start: "AAT", offset: 40}, "prod2" => {start: "TGG", offset: 10} }})
        .keys)
        .to include("AAT", "TGG") #: {start: "TGG", offset: 10} }) 
    end

    it "checks all product starts are the same length" do
      expect(MutationsCLI.construct_products({ products: {"prod1" => {start: "AAT", offset: 40}, "prod2" => {start: "TGG", offset: 10} }})[:start_length])
        .to eq 3
    end
  end

  describe "Calling mutations" do
    it "Calls mutations from a SAM alignment given a config hash" do
      sam = Bio::DB::Alignment.new("I2M5K:00271:01406	0	5	112839854	70	55M12D162M7S	*	0	0	CAGTGATCTTCCAGATAGCCCTGGACAAACCATGCCACCAAGCAGAAGTAAAACACCTCAAACAGCTCAAACCAAGCGAGAAGTACCTAAAAATAAAGCACCTACTGCTGAAAAGAGAGAGAGTGGACCTAAGCAAGCTGCAGTAAATGCTGCAGTTCAGAGGGTCCAGGTTCTTCCAGATGCTGATACTTTATTACATTTTGCCACGGAAAGTACTGCTGAGG	ACHECCC@???ACCCCCCDC>CC@CDCC>C>>?CC=>?ACADCCCCACCCCC:AAA<=CCC??<>CDCDE?C@C>=;CC=>>>@@@>:::::+:::/:@@@<>>?CCCD=>>=7:D???AAAAAB8;;>??;?@@=?;;;???@@@:@B;;;;GBB?BBAAAA9??=@@;?<?A>B?C@@@@@CBBB?;;;4;C?BB;;;:1:B>BBB=AA;A@@???A>?::2	PG:Z:novoalign	AS:i:194	UQ:i:194	NM:i:12	MD:Z:55^CCTCCACCACCT162")
      config = { 
          start: "CAG",
          offset: 30,
          length: 50,
          translation_start: 1
        }
      mutation_array = MutationsCLI.call_mutations_given_product(sam,config)
      expect(mutation_array.to_hgvs).to eq "5:g.112839909_112839920delCCTCCACCACCT"
    end
  end


end
