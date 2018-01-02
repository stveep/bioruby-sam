require 'helper'
class SAMTest < Test::Unit::TestCase
	def test_split_types
		sam = Bio::DB::Alignment.new("DKNQZ:00025:00303	0	5	112767204	37	60M1D7M2I6M	*	0	0	GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA	CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:3	XO:i:1	XG:i:1	MD:Z:60^G13")
		assert(sam.qname.is_a? String)
		assert(sam.flag.is_a? Integer) # DB::Alignment only supports integer FLAG
		assert(sam.rname.is_a? String)
		assert(sam.pos.is_a? Integer)
		assert(sam.mapq.is_a? Integer)
		assert(sam.cigar.is_a? String)
		assert(sam.mrnm.is_a?(String) || sam.mrnm.nil?)
		assert(sam.mpos.is_a?(Integer) || sam.mpos.nil?)
		assert(sam.isize.is_a?(Integer) || sam.isize.nil?)
		assert(sam.seq.is_a? String)
		assert(sam.qual.is_a? String)
		assert(sam.tags.is_a? Hash)

	end
	def test_split
		sam = Bio::DB::Alignment.new("DKNQZ:00025:00303	0	5	112767204	37	60M1D7M2I6M	*	0	0	GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA	CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:3	XO:i:1	XG:i:1	MD:Z:60^G13")
		assert_equal(sam.qname,"DKNQZ:00025:00303","ID not as expected")
		assert_equal(sam.flag,0,"Flag not as expected")
		assert_equal(sam.rname,"5","Chr not as expected")
		assert_equal(sam.pos,112767204,"Position not as expected")
		assert_equal(sam.mapq,37,"Quality not as expected")
		assert_equal("60M1D7M2I6M", sam.cigar)
		assert_equal("*",sam.mrnm)
		assert_equal(0,sam.mpos)
		assert_equal(0,sam.isize)
		assert_equal(sam.seq,"GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA")
		assert_equal(sam.qual,"CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112","Base quality string not as expected")
		assert_equal("60^G13",sam.tags["MD"].value)

	end

	def test_aliases
		sam = Bio::DB::Alignment.new("DKNQZ:00025:00303	0	5	1	37	60M1D7M2I6M	*	0	0	GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA	CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:3	XO:i:1	XG:i:1	MD:Z:60^G13")
		assert_equal(sam.chr,sam.rname)
		assert_equal(sam.tags,sam.opt)
	end

	def test_insertion_mutation
		insertion = Bio::DB::Alignment.new("I2M5K:00253:00406	0	5	112839854	70	63M2I138M1D27M7S	*	0	0	CAGTGATCTTCCAGATAGCCCTGGACAAACCATGCCACCAAGCAGAAGTAAAACACCTCCACCATACCTCCTCAAACAGCTCAAACCAAGCGAGAAGTACCTAAAAATAAAGCACCTACTGCTGAAAAGAGAGAGAGTGGACCTAAGCAAGCTGCAGTAAATGCTGCAGTTCAGAGGGTCCAGGTTCTTCCAGATGCTGATACTTATTACATTTTGCCACGGAAAGTACTGCTGAGG	@CDDDCCCCACACCCCCCCC?CCACCCC>A6;;;;7;;6;6;BC;;6;;;;;.;;>ADDA??;;;;;?CCACCCD>C??@CCCC>C@C;>?CCCC@C=::@:::::+:::/:CCC?>>>>CCCCDDD9CCCC@AB????=AB>??;?BB>@@@AA???CC<@@?????BB>??;;;B<BC;??8;6:A=@=@BBB;;;?<77//*08*088888*8=9=?B7;;4;??????????<	PG:Z:novoalign	AS:i:183	UQ:i:183	NM:i:3	MD:Z:201^T27")
		assert_equal 112839916, insertion.mutations[0].position
		assert_nil insertion.mutations[0].reference
		assert_equal "AT", insertion.mutations[0].mutant
		assert_equal :insertion, insertion.mutations[0].type
		assert_equal 112839916, insertion.mutations(60,10)[0].position
	end
  def test_ins_with_offset
		insertion = Bio::DB::Alignment.new("I2M5K:00253:00406	0	5	112839854	70	63M2I138M1D27M7S	*	0	0	CAGTGATCTTCCAGATAGCCCTGGACAAACCATGCCACCAAGCAGAAGTAAAACACCTCCACCATACCTCCTCAAACAGCTCAAACCAAGCGAGAAGTACCTAAAAATAAAGCACCTACTGCTGAAAAGAGAGAGAGTGGACCTAAGCAAGCTGCAGTAAATGCTGCAGTTCAGAGGGTCCAGGTTCTTCCAGATGCTGATACTTATTACATTTTGCCACGGAAAGTACTGCTGAGG	@CDDDCCCCACACCCCCCCC?CCACCCC>A6;;;;7;;6;6;BC;;6;;;;;.;;>ADDA??;;;;;?CCACCCD>C??@CCCC>C@C;>?CCCC@C=::@:::::+:::/:CCC?>>>>CCCCDDD9CCCC@AB????=AB>??;?BB>@@@AA???CC<@@?????BB>??;;;B<BC;??8;6:A=@=@BBB;;;?<77//*08*088888*8=9=?B7;;4;??????????<	PG:Z:novoalign	AS:i:183	UQ:i:183	NM:i:3	MD:Z:201^T27")
		assert_equal 112839916, insertion.mutations(60,10)[0].position
	end
	def test_sorting
		insertion_and_deletion = Bio::DB::Alignment.new("I2M5K:00253:00406	0	5	1	70	63M2I138M1D27M7S	*	0	0	CAGTGATCTTCCAGATAGCCCTGGACAAACCATGCCACCAAGCAGAAGTAAAACACCTCCACCATACCTCCTCAAACAGCTCAAACCAAGCGAGAAGTACCTAAAAATAAAGCACCTACTGCTGAAAAGAGAGAGAGTGGACCTAAGCAAGCTGCAGTAAATGCTGCAGTTCAGAGGGTCCAGGTTCTTCCAGATGCTGATACTTATTACATTTTGCCACGGAAAGTACTGCTGAGG	@CDDDCCCCACACCCCCCCC?CCACCCC>A6;;;;7;;6;6;BC;;6;;;;;.;;>ADDA??;;;;;?CCACCCD>C??@CCCC>C@C;>?CCCC@C=::@:::::+:::/:CCC?>>>>CCCCDDD9CCCC@AB????=AB>??;?BB>@@@AA???CC<@@?????BB>??;;;B<BC;??8;6:A=@=@BBB;;;?<77//*08*088888*8=9=?B7;;4;??????????<	PG:Z:novoalign	AS:i:183	UQ:i:183	NM:i:3	MD:Z:201^T27")
		assert_equal :insertion, insertion_and_deletion.mutations[0].type
		assert_equal :deletion, insertion_and_deletion.mutations[1].type
	end

	def test_deletion_mutation
		deletion = Bio::DB::Alignment.new("I2M5K:00271:01406	0	5	112839854	70	55M12D162M7S	*	0	0	CAGTGATCTTCCAGATAGCCCTGGACAAACCATGCCACCAAGCAGAAGTAAAACACCTCAAACAGCTCAAACCAAGCGAGAAGTACCTAAAAATAAAGCACCTACTGCTGAAAAGAGAGAGAGTGGACCTAAGCAAGCTGCAGTAAATGCTGCAGTTCAGAGGGTCCAGGTTCTTCCAGATGCTGATACTTTATTACATTTTGCCACGGAAAGTACTGCTGAGG	ACHECCC@???ACCCCCCDC>CC@CDCC>C>>?CC=>?ACADCCCCACCCCC:AAA<=CCC??<>CDCDE?C@C>=;CC=>>>@@@>:::::+:::/:@@@<>>?CCCD=>>=7:D???AAAAAB8;;>??;?@@=?;;;???@@@:@B;;;;GBB?BBAAAA9??=@@;?<?A>B?C@@@@@CBBB?;;;4;C?BB;;;:1:B>BBB=AA;A@@???A>?::2	PG:Z:novoalign	AS:i:194	UQ:i:194	NM:i:12	MD:Z:55^CCTCCACCACCT162")
		assert_equal deletion.mutations[0].type, :deletion
		assert_equal deletion.mutations[0].reference, "CCTCCACCACCT"
		assert_equal deletion.mutations[0].mutant, nil
		assert_equal deletion.mutations(50,40)[0].position, 112839909
		assert_equal "5", deletion.mutations[0].seqname
	end
	def test_substitution_mutation
		substitution = Bio::DB::Alignment.new("OR1FQ:00462:02257	0	ENST00000366794	936	70	193M43S	*	0	0	ACTCATCTTCAACAAGCAGCAAGTGCCTTCTGGGGAGTCGGCGATCTTGGACCGAGTAGCCGATGGCATGGTGTTCGGTGCCCTCCTTCCCTGCGAGGAATGCTCGGGTCAGCTGGTCTTCAAGAGCGATGCCTATTACTGCACTGGGGACGTCACTGCCTGGACCAAGTGTATGGTCAAGACACAGACACCCTCCACAGCCTCGGCTCCTGCTGCTGTGAACTCCTCTGCTGAGG	@BCACC@@@@DAFFADCCCCDEID@@?@ACCCDD:@??;<8<1..=>=<1111@@CD??@@CC@C@CFDCACCCADDABCCD?DD@CACD?CC??>C6;6;>>???E?C@??CCDACCC@CD@CCC><?A>>7;<<7<<<??;BBBBB/;;;;;BBBCC@CCACC=@7;;;;;;;6;@@7?@CCCCCC111<6<<+00>>>=CCD;??C@CCCC?????CCAC????CCECDDCDB	PG:Z:novoalign	AS:i:328	UQ:i:328	NM:i:1	MD:Z:60T132")
		assert_equal substitution.mutations[0].type, :substitution
		assert_equal substitution.mutations[0].reference, "T"
		assert_equal substitution.mutations[0].mutant, "C"
		assert_equal 996, substitution.mutations[0].position
		assert_equal "ENST00000366794", substitution.mutations[0].seqname
	end
	def test_substitution_mutation_with_translation_pos
		substitution = Bio::DB::Alignment.new("OR1FQ:00462:02257	0	ENST00000366794	936	70	193M43S	*	0	0	ACTCATCTTCAACAAGCAGCAAGTGCCTTCTGGGGAGTCGGCGATCTTGGACCGAGTAGCCGATGGCATGGTGTTCGGTGCCCTCCTTCCCTGCGAGGAATGCTCGGGTCAGCTGGTCTTCAAGAGCGATGCCTATTACTGCACTGGGGACGTCACTGCCTGGACCAAGTGTATGGTCAAGACACAGACACCCTCCACAGCCTCGGCTCCTGCTGCTGTGAACTCCTCTGCTGAGG	@BCACC@@@@DAFFADCCCCDEID@@?@ACCCDD:@??;<8<1..=>=<1111@@CD??@@CC@C@CFDCACCCADDABCCD?DD@CACD?CC??>C6;6;>>???E?C@??CCDACCC@CD@CCC><?A>>7;<<7<<<??;BBBBB/;;;;;BBBCC@CCACC=@7;;;;;;;6;@@7?@CCCCCC111<6<<+00>>>=CCD;??C@CCCC?????CCAC????CCECDDCDB	PG:Z:novoalign	AS:i:328	UQ:i:328	NM:i:1	MD:Z:60T132")
		# offset, length, reference start, translation start
		assert_equal 852, substitution.mutations(55,20,145)[0].position
		assert_equal "ENST00000366794", substitution.mutations[0].seqname
	end

	def test_cdna_mutation
		skip_unless_ensembl_available do
			parp1	=	Bio::DB::Alignment.new("OR1FQ:00028:00030	0	ENST00000366794	342	70	66M1D88M6D70M8S	*	0	0	CCCTGACGTTGAGGTGGATGGGTTCTCTGAGCTTCGGTGGGATGATCAGCAGAAAGTCAAGAAGACGCGGAAGCTGGAGGAGTGACAGGCAAAGGCCAGGATGGAATTGGTAGCAAGGCAGAGAAGACTCTGGGTGACTTTGCAGCAGAGTATGCCAACAGAAGTACGTGCAAGGGGTGTATGGAGAAGATAGAAAAGGGCCAGGTGCGCCTGTCCAAGAAGATGGCTGAGG	;;1;;;;6606660;B?A<<<1?ACCDC?@;;;A<;7;<<16B==BDB@@@;;;1;@@@:;/*/;/0--)-)-C660>B@=?@D?;;;7;;;1;7;@;;7;;;64.4.4.454;;6;=@CFDCC@?>;@A;;>:;CACCC>CCCCCCCCCCCCCC@C>@@CCACCCCECC@@6:::.::::>D>?>CACEC?>1<<(00*0*0/6777?A??C??6;?;;6;@::;?CDCD>	PG:Z:novoalign	AS:i:234	UQ:i:234	NM:i:8	MD:Z:45C20^A88^CCAAGT70")
			response = parp1.mutations(140,40,145).vep("human","c")
			assert_equal "ENST00000366794:c.353_358delCCAAGT", parp1.mutations(140,40,145).first.to_hgvs("c")
			assert_equal "CCAAGT/-", response.first["allele_string"]
			assert_equal 226392243, response.first["start"]
  	end
	end

	def test_nil_if_no_mutation
		no_mut = Bio::DB::Alignment.new("DKNQZ:00025:00303	0	5	1	37	75M	*	0	0	GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA	CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:3	XO:i:1	XG:i:1	MD:Z:75")
	  assert_nil no_mut.mutations
	end

	def test_query_no_mutation
		no_mut = Bio::DB::Alignment.new("DKNQZ:00025:00303	0	5	1	37	75M	*	0	0	GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA	CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:3	XO:i:1	XG:i:1	MD:Z:75")
		assert_equal "TTTCCCTGGA", no_mut.query(8,10)
	end
	def test_query_deletion
		del_mut = Bio::DB::Alignment.new("DKNQZ:00025:00303	0	5	100000	37	9M1D65M	*	0	0	GCAGTAATTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA	CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:3	XO:i:1	XG:i:1	MD:Z:9^T65")
		assert_equal "TT-CCCTGGA", del_mut.query(8,10)
	end
	def test_query_insertion
		ins_mut = Bio::DB::Alignment.new("DKNQZ:00025:00303	0	5	100000	37	10M2I65M	*	0	0	GCAGTAATTTGGCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA	CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:3	XO:i:1	XG:i:1	MD:Z:75")
		assert_equal "TTT_gg_CCCTGGA", ins_mut.query(8,10)
	end
	def test_query_substitution
		sub_mut = Bio::DB::Alignment.new("DKNQZ:00025:00303	0	5	100000	37	75M	*	0	0	GCAGTAATTTCGCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA	CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:3	XO:i:1	XG:i:1	MD:Z:11C63")
		assert_equal "TTTCgCTGGA", sub_mut.query(8,10)
	end

	def test_query_complex
	complex	=	Bio::DB::Alignment.new	"MLF7W:00161:00415	0	1	226064342	70	121M30I32M1D52M7S	*	0	0	TCTCTGTACCATGGCTCGTACAAAGCAGACTGCCCGCAAATCGACCGGTGGTAAAGCACCCAGGAAGCAACTGGCTACAAAAGCCGCTCGCAAGAGTGCGCCCTCTACTGGAGTGGTGAAGAAACCTCATCGTTAAAAAAACCTCATCGTTAAACCTCATCGTTACAGGTATTAAAAAACAGGAAAAAATGGGACAAAGTCTCTCTTGTATGTATCCATATAATTTAACAAAAAGAATGGAT	CC???ECCCACCCACBBCD>?@D=@@DCEEDDCC=CCDD;??;;;=@=;;0;66,65;5;1;@?C@CCCACCC=CCCCCCB.;;6;;;C;=:@:9:=2::>>52:BBCBCAD?>=ADCC?==:/:><<<=8884:CCCCC(B288::::C3=@8C?@:8888/*//0-'--.'-:;:;,><=288888)6BE?DCCC?DC@@>@;;7;00000<<DB@@@DC@CE?CA00000(0/*/;+;0	PG:Z:novoalign	AS:i:379	UQ:i:379	NM:i:33	MD:Z:113G39^A35C16	YH:m:1:g.[226064455G>T;226064462_226064463insAAACCTCATCGTTAAAAAAACCTCATCGTT;226064495delA]"
		# offset 100, length 70
		# MD:Z:113G39^A35C16      YH:m:1:g.[226064455G>T;226064462_226064463insAAACCTCATCGTTAAAAAAACCTCATCGTT;226064495delA]
		assert_equal "GCCCTCTACTGGAGtGGTGAAG_aaacctcatcgttaaaaaaacctcatcgtt_AAACCTCATCGTTACAGGTATTAAAAAACAGG-", complex.query(100,70)
	end

	def test_query_with_preceding
	pre = Bio::DB::Alignment.new	"MLF7W:00741:00261	0	1	226064342	70	32M1I127M2D33M1D34M1I14M36S	*	0	0	TCTCTGTACCATGGCTCGTACAAAGCAGACTGCCCCGCAAATCGACCGGTGGTAAAGCACCCAGGAAGCAACTGGCTACAAAAGCCGCTCGCAAGAGTGCGCCCTCTACTGGAGTGGTGAAGAAACCTCATCGTTACAGGTATTAAAAAACAGGAAAAAATGGACAAAGTCTCTCTTGTATGTATCCATATAATTAACAAAAAGATGGATAACAGGAAAACTTTTTGTCTTTAGAGAAACTTCATTTTGAACACTTAAACTTTACTGCTTAGCTGAGG	CC?>=CB@B>?@C@CCCDBB<66,<AA;;;;;;;;,6=CC,>>>>>@C@BC@;;;1;>;>C????@???DA???ACCCCCDC8@@ACCCFC??>?CCCCCCCC<CCCCCCE@CCCCACDCACCE:>=;:80:00*////)/??288888)8005:::::(0.).-61)133::BBC;????<<<882888?D20*/*/////(..66-.-------999-9;;99,0...8-:6:/88,8<7>?ABB1;:9888<288,888,88=88/*//0008:*	PG:Z:novoalign	AS:i:534	UQ:i:534	NM:i:8	MD:Z:113G45^AT0G27C4^T48	YH:m:1:g.[226064455G>T;226064501_226064502delAT;226064503G>T]"
	# MD:Z:113G45^AT0G27C4^T48        YH:m:1:g.[226064455G>T;226064501_226064502delAT;226064503G>T]
	assert_equal "GCCCTCTACTGGAGtGGTGAAGAAACCTCATCGTTACAGGTATTAAAAAACAGGAAAAAA--tGGACAAA", pre.query(100,70)
	end

	def test_query_multi_sub
		rpts =  Bio::DB::Alignment.new("OR1FQ:00079:00560\t0\tENST00000366794\t936\t70\t361M4S\t*\t0\t0\tACTCATCTTCAACAAGCAGCAAGTGCCTTCTGGGGAGTCGGCGATCTTGGACCGAGTAGCCGATGGCATGGTGTTCGGTGCCCTCCTTCCCTGCGAGGAATGCTCGGGTCAGCTGGTCTTCAAGAGCGATGCCTATTACTGCACTGGGGACGTCACTGCCTGGACCAAGTGTATGGTCAAGACACAGACACCCCAACCGAAGGAGTGGGTAACCCCAAAGGAATTCCGAGAAATCTCTTACCTCAAGAAATTGAAGGTTAAAAAACAGGACCGTATATTCCCCCCAGAAACCAGCGCCTCCGTGGCGGCCACGCCTCCGCCCTCCACAGCCTCGGCTCCTGCTGCTGTGAACTCCTCTGCTGAGG\tCCDCCCBCACCAECADCDCCDACCCDAFAC?;;;/;;>BCAB;;;;;7;?BB>@>CC;;>7;B@@=@;<;7?CCAEEADCCC?@@?@ADD?CCCCCC6>>???CCEE>C???CCCBCCD???@CCCD?><;<7<A@ACHDDCCCCCCC;C??>B>;;CB>@B@CC@CADFC???FAD??>????CCCCC;;;;.;6;6;@6;6>>00<600*///(-5/537=DAC@BC;;>:;;;AC>?@AECDACCD>DADD@B@BBDDDDD1CCCACC6;;;;;AACCDCD1C>??<?ACFDDDACC?;;A><<7<7<;;;ADDACCD?DF@@@DCC>CCCADDHACDDEDC@@@CB>A;;7;;;;CDFCCA\tPG:Z:novoalign\tAS:i:160\tUQ:i:160\tNM:i:4\tMD:Z:60T132A1C1G163")

		assert_equal "ACTCATCTTCAACAAGCAGCAAGTGCCTTCTGGGGAGTCGGCGATCTTGGACCGAGTAGCcGATGGCATGGTGTTCGGTGCCCTCCTTCCCTGCGAGGAATGCTCGGGTCAGCTGGTCTTCAAGAGCGATGCCTATTACTGCACTGGGGACGTCACTGCCTGGACCAAGTGTATGGTCAAGACACAGACACCCcAaCcGAAGGAGTGGGTAACCCCAAAG",
		rpts.query(1,220)
	end


	def test_complex_mutation
		sub_del_ins = Bio::DB::Alignment.new "OR1FQ:00021:00043	0	ENST00000366794	936	70	128M1D16M1I38M	*	0	0	ACTCATCTTCAACAAGCAGCAAGTGCCTTCTGGGGAGTCGGCGATCTTGGACCGAGTAGCCGATGGCATGGTGTTCGGTGCCCTCCTTCCCTGCGAGGAATGCTCGGGTCAGCTGGTCTTCAAGAGCGTGCCTATTACTGCACTGGGGGACGTCACTGCCTGGACCAAGTGTATGGTCAAGAC	CCCCCCCDB@?=<B7;;<<<A7?@FCDDABBCCD;C???BAB?@@?CAC@CC>??C@;;;7;;;/8/000+.;;8/@7<;;;1;;7;7;;1>BBCCDAD;???;;;C:;;C;;;?@@BC=;;7;;<00...)-55357DC<<6;;;;;,66;;;;;;6606606<:@@;5;44--)--.)---	PG:Z:novoalign	AS:i:121	UQ:i:121	NM:i:3	MD:Z:60T67^A54"
		# MD:Z:60T67^A54
		# Old "manual" method:
		assert_equal "996T>C;1064delA;1080_1081insG", sub_del_ins.mutations.map{|m| m.to_hgvs}.join(";")
	end

	def test_really_complex_mutation
		omg = Bio::DB::Alignment.new "MLF7W:01389:01808	0	ENST00000366794	957	70	153M45D46M1I83M1I13M4S	*	0	0	AGTGCCTTCTGGGGAGTCGGCGATCTTGGACCGAGTAGCCGATGGCATGGTGTTCGGTGCCCTCCTTCCCTGCGAGGAATGCTCGGGTCAGCTGGTCTTCAAGAGCGATGCCTATTACTGCACTGGGGACGTCACTGCCTGGACCAAGTGTATGGAATTCCGAGAAATCTCTTACCTCAAGAAATTGAAGGTTAAAAAAGCAGGACCGTATATTCCCCCCAGAAACCAGCGCCTCCGTGGCGGCCACGCCTCCGCCCTCCACAGCCTCGGCTCCTGCTGCTGTCGAACTCCTCTGCTGAGG	0666606066666,66@//*//----)-)--)666@6666B@@<666660---)--)----)660606C<65>5;;6;4;BCCCCD>CCC??>B7;;C7;;@AAAA<<<<<7<<<7<;DA;;;;---%-------6BB6??=;D<@7;;;;;;A6;7;7;7;;<6D:6<B<<7;;7;;;7;;;1<><;7;@C7;;;;;+;BCCA?><<;BB;;7;;;;;+;;;;160666.-)-.)--6)-4.;6;;AB>;;7;;;1;C7;;5;;0666066<0555<?6;;/.-.).8==CCCCCA;:/6	PG:Z:novoalign	AS:i:504	UQ:i:504	NM:i:48	MD:Z:39T113^GGTCAAGACACAGACACCCAACCGGAAGGAGTGGGTAACCCCAAA142"
		# MDZ MD:Z:39T113^GGTCAAGACACAGACACCCAACCGGAAGGAGTGGGTAACCCCAAA142
		# translation start 145
		assert_equal "ENST00000366794:c.[852T>C;966_1010delGGTCAAGACACAGACACCCAACCGGAAGGAGTGGGTAACCCCAAA;1056_1057insG;1139_1140insC]", omg.mutations(1,nil,145).to_hgvs

		#This one (wrong mutant base retrieved for substition) seems to be a problem with novoalign? Mutation is CT>TC but aligned as a single T substition
		#wtf = Bio::DB::Alignment.new "MLF7W:00090:00609\t0\tENST00000366794\t957\t70\t25M1I314M4S\t*\t0\t0\tAGTGCCTTCTGGGGAGTCGGCGATCTTTGGACCGAGTAGTCGATGGCATGGTGTTCGGTGCCCTCCTTCCCTGCGAGGAATGCTCGGGTCAGCTGGTCTTCAAGAGCGATGCCTATTACTGCACTGGGGACGTCACTGCCTGGACCAAGTGTATGGTCAAGACACAGACACCCAACCGGAAGGAGTGGGTAACCCCAAAGGAATTCCGAGAAATCTCTTACCTCAAGAAATTGAAGGTTAAAAAGCAGGACCGTATATTCCCCCCAGAAACCAGCGCCTCCGTGGCGGCCACGCCTCCGCCCTCCACAGCCTCGGCTCCTGCTGCTGTGAACTCCTCTGCGAGG\t?AA<<7;7>>DDC;@@@DCAB;;;;;;1;7;@?;;;;;;;6;;;@@C;<<7<<;7<?7<CBB:<;7;7;;1;;B;;;7;=@BCCCCC1>>;;>BC@CCDACDACC;;<<<<<7?FEADCDGDDA;?;;/;<<;;CCB@CA@@>?DAC???DCCFDAD??>?CCCDDCD@??C9@<B;>=?@C@CFCCD>CC6;;;/;C;CAD>?>?ACCCDC?DDCCCACCACCCA@@@:>6;E>A@C<;;;;-;CCCADD<;;C>??ACGGGE+;<<<1;>ACCDDADC7<;B7/0+0+////8;;7?AA9<<66<6;<66CC?C?><;<<;;;;>;ACA>??A>>;;>>;;7\tPG:Z:novoalign\tAS:i:164\tUQ:i:164\tNM:i:3\tMD:Z:39T203A95"
		# MD:Z:39T203A95
		#assert_equal "ENST00000366794:c.[981_982insT;996T>T;1200A>A]", wtf.mutations.to_hgvs

	end

	def test_add_tag
		sub_del_ins = Bio::DB::Alignment.new "OR1FQ:00021:00043	0	ENST00000366794	936	70	128M1D16M1I38M	*	0	0	ACTCATCTTCAACAAGCAGCAAGTGCCTTCTGGGGAGTCGGCGATCTTGGACCGAGTAGCCGATGGCATGGTGTTCGGTGCCCTCCTTCCCTGCGAGGAATGCTCGGGTCAGCTGGTCTTCAAGAGCGTGCCTATTACTGCACTGGGGGACGTCACTGCCTGGACCAAGTGTATGGTCAAGAC	CCCCCCCDB@?=<B7;;<<<A7?@FCDDABBCCD;C???BAB?@@?CAC@CC>??C@;;;7;;;/8/000+.;;8/@7<;;;1;;7;7;;1>BBCCDAD;???;;;C:;;C;;;?@@BC=;;7;;<00...)-55357DC<<6;;;;;,66;;;;;;6606606<:@@;5;44--)--.)---	PG:Z:novoalign	AS:i:121	UQ:i:121	NM:i:3	MD:Z:60T67^A54"
		sub_del_ins.add_tag!("TE:s:TTAG")
		assert_equal "TTAG", sub_del_ins.tags["TE"].value

		tag_obj = Bio::DB::Tag.new
		tag_obj.set("TA:g:Object")
		sub_del_ins.add_tag!(tag_obj)
		assert_equal "Object", sub_del_ins.tags["TA"].value

		 assert_match "TA:g:Object", sub_del_ins.sam_string
	end

end
