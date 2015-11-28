require 'helper'
class SAMTest < Test::Unit::TestCase
	def test_split_types
		sam = Bio::Alignment::SAM.new("DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U NM:i:3 X0:i:1 X1:i:0 XM:i:3 XO:i:1 XG:i:1 MD:Z:60^G13")
		assert(sam.qname.is_a? String)
		assert(sam.flag.is_a? String)
		assert(sam.rname.is_a? String)
		assert(sam.pos.is_a? Integer)
		assert(sam.mapq.is_a? Integer)
		assert(sam.cigar.instance_of? Bio::Alignment::CIGAR)
		assert(sam.mrnm.is_a?(String) || sam.mrnm.nil?)
		assert(sam.mpos.is_a?(Integer) || sam.mpos.nil?)
		assert(sam.isize.is_a?(Integer) || sam.isize.nil?)
		assert(sam.seq.instance_of? Bio::Sequence::NA)
		assert(sam.qual.is_a? String)
		assert(sam.tags.is_a? Hash)

	end
	def test_split
		sam = Bio::Alignment::SAM.new("DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U  NM:i:3  X0:i:1  X1:i:0  XM:i:3  XO:i:1  XG:i:1  MD:Z:60^G13")
		assert_equal(sam.qname,"DKNQZ:00025:00303","ID not as expected")
		assert_equal(sam.flag,"0","Flag not as expected")
		assert_equal(sam.rname,"5","Chr not as expected")
		assert_equal(sam.pos,112767204,"Position not as expected")
		assert_equal(sam.mapq,37,"Quality not as expected")
		assert_equal(sam.cigar.pairs.length,5,"CIGAR not as expected")
		assert_equal(sam.mrnm,nil,"Mate chr not as expected [* -> nil]")
		assert_equal(sam.mpos,nil,"Mate pos not as expected [0 -> nil]")
		assert_equal(sam.isize,nil,"Insertion size not as expected [0 -> nil]")
		assert_equal(sam.seq,Bio::Sequence::NA.new("GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA"),"Sequence not as expected")
		assert_equal(sam.qual,"CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112","Base quality string not as expected")
 assert_equal(sam.tags,{:unsplit_tags => "XT:A:U NM:i:3 X0:i:1 X1:i:0 XM:i:3 XO:i:1 XG:i:1 MD:Z:60^G13", "XT:A" => "U", "NM:i" => "3", "X0:i" => "1", "X1:i" => "0", "XM:i" => "3", "XO:i" => "1","XG:i" => "1", "MD:Z" =>"60^G13"})

	end

	def test_aliases
		sam = Bio::Alignment::SAM.new("DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U  NM:i:3  X0:i:1  X1:i:0  XM:i:3  XO:i:1  XG:i:1  MD:Z:60^G13")
		assert_equal(sam.chr,sam.rname)
		assert_equal(sam.tags,sam.opt)
	end

	def test_insertion_mutation
		insertion = Bio::Alignment::SAM.new("I2M5K:00253:00406       0       5       112839854       70      63M2I138M1D27M7S        *       0       0       CAGTGATCTTCCAGATAGCCCTGGACAAACCATGCCACCAAGCAGAAGTAAAACACCTCCACCATACCTCCTCAAACAGCTCAAACCAAGCGAGAAGTACCTAAAAATAAAGCACCTACTGCTGAAAAGAGAGAGAGTGGACCTAAGCAAGCTGCAGTAAATGCTGCAGTTCAGAGGGTCCAGGTTCTTCCAGATGCTGATACTTATTACATTTTGCCACGGAAAGTACTGCTGAGG   @CDDDCCCCACACCCCCCCC?CCACCCC>A6;;;;7;;6;6;BC;;6;;;;;.;;>ADDA??;;;;;?CCACCCD>C??@CCCC>C@C;>?CCCC@C=::@:::::+:::/:CCC?>>>>CCCCDDD9CCCC@AB????=AB>??;?BB>@@@AA???CC<@@?????BB>??;;;B<BC;??8;6:A=@=@BBB;;;?<77//*08*088888*8=9=?B7;;4;??????????<   PG:Z:novoalign  AS:i:183        UQ:i:183        NM:i:3  MD:Z:201^T27")
		assert_equal :insertion, insertion.mutations[0].type
		assert_nil insertion.mutations[0].reference
		assert_equal "AT", insertion.mutations[0].mutant
		assert_equal 63, insertion.mutations[0].position
	end
	def test_deletion_mutation
  	deletion = Bio::Alignment::SAM.new("I2M5K:00271:01406       0       5       112839854       70      55M12D162M7S    *       0       0       CAGTGATCTTCCAGATAGCCCTGGACAAACCATGCCACCAAGCAGAAGTAAAACACCTCAAACAGCTCAAACCAAGCGAGAAGTACCTAAAAATAAAGCACCTACTGCTGAAAAGAGAGAGAGTGGACCTAAGCAAGCTGCAGTAAATGCTGCAGTTCAGAGGGTCCAGGTTCTTCCAGATGCTGATACTTTATTACATTTTGCCACGGAAAGTACTGCTGAGG        ACHECCC@???ACCCCCCDC>CC@CDCC>C>>?CC=>?ACADCCCCACCCCC:AAA<=CCC??<>CDCDE?C@C>=;CC=>>>@@@>:::::+:::/:@@@<>>?CCCD=>>=7:D???AAAAAB8;;>??;?@@=?;;;???@@@:@B;;;;GBB?BBAAAA9??=@@;?<?A>B?C@@@@@CBBB?;;;4;C?BB;;;:1:B>BBB=AA;A@@???A>?::2        PG:Z:novoalign  AS:i:194        UQ:i:194        NM:i:12 MD:Z:55^CCTCCACCACCT162")
		assert_equal deletion.mutations[0].type, :deletion
		assert_equal deletion.mutations[0].reference, "CCTCCACCACCT"
		assert_equal deletion.mutations[0].mutant, nil
		assert_equal deletion.mutations[0].position, 56
  end
	def test_substitution_mutation
		substitution = Bio::Alignment::SAM.new("OR1FQ:00462:02257       0       ENST00000366794 936     70      193M43S *       0       0       ACTCATCTTCAACAAGCAGCAAGTGCCTTCTGGGGAGTCGGCGATCTTGGACCGAGTAGCCGATGGCATGGTGTTCGGTGCCCTCCTTCCCTGCGAGGAATGCTCGGGTCAGCTGGTCTTCAAGAGCGATGCCTATTACTGCACTGGGGACGTCACTGCCTGGACCAAGTGTATGGTCAAGACACAGACACCCTCCACAGCCTCGGCTCCTGCTGCTGTGAACTCCTCTGCTGAGG    @BCACC@@@@DAFFADCCCCDEID@@?@ACCCDD:@??;<8<1..=>=<1111@@CD??@@CC@C@CFDCACCCADDABCCD?DD@CACD?CC??>C6;6;>>???E?C@??CCDACCC@CD@CCC><?A>>7;<<7<<<??;BBBBB/;;;;;BBBCC@CCACC=@7;;;;;;;6;@@7?@CCCCCC111<6<<+00>>>=CCD;??C@CCCC?????CCAC????CCECDDCDB    PG:Z:novoalign  AS:i:328        UQ:i:328        NM:i:1  MD:Z:60T132")
		assert_equal substitution.mutations[0].type, :substitution
		assert_equal substitution.mutations[0].reference, "CCTCCACCACCT"
		assert_equal substitution.mutations[0].mutant, nil
		assert_equal substitution.mutations[0].position, 56
  end
end
