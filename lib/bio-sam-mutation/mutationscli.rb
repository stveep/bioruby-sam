module MutationsCLI
  @default_species = "human"
  @tag_to_add = "YH:m"
  @comment_char = "#"
  @barcode_regexp = /IonXpress_(\w+)\.R/
  class << self
    attr_accessor :default_species, :tag_to_add, :comment_char, :barcode_regexp
  end

  def self.tag (sam, config, species=MutationsCLI.default_species)
    # With multiple files, could end up with duplicate headers...
    if sam.match(/^@/)
      if config[:single_product]
        config[:outfile].puts sam
        return
      else
        # If we have multiple products, output sam headers to each output file
        config.each do |key, config_hash|
          next if [:start_length].include? key
          config_hash[:outfile].puts sam
        end
        return
      end
    end
    sam = Bio::DB::Alignment.new(sam)
    # DRY up:
    unless config[:single_product]
      first_bases = sam.seq[0..config[:start_length]-1].upcase
      config[first_bases] ? config = config[first_bases] : return
    end
    # Stop if a search file is specified, and this isn't it.
    return if MutationsCLI.not_included_file?(config, ARGF)
    mut_array = sam.mutations(config[:offset],config[:length],config[:translation_start])
    if mut_array
      new_tag = Bio::DB::Tag.new(MutationsCLI.tag_to_add + ":" + mut_array.to_hgvs)
      config[:outfile].puts sam.add_tag(new_tag)
    end
  end

  # Returns a MutationArray
  def self.call_mutations sam, config
    unless config[:single_product]
      first_bases = sam.seq[0..config[:start_length]-1].upcase
      config[first_bases] ? config = config[first_bases] : return
    end
    MutationsCLI.call_mutations_given_product sam, config
  end

  def self.call_mutations_given_product sam, config
    # Stop if a search file is specified, and this isn't it.
    return if MutationsCLI.not_included_file?(config, ARGF)
    sam.mutations(config[:offset],config[:length],config[:translation_start])
  end

  def self.not_included_file? config, input
    config[:file] && input.filename !~ config[:file]
  end

  def self.set_defaults config_hash
    config_hash[:start] ||= "." # i.e. regexp will match anything
    config_hash[:offset] ||= 1
    config_hash[:length] ||= 100
    config_hash[:translation_start] ||= 1
    config_hash[:file] = Regexp.new(config_hash[:file]) if config_hash[:file]
    config_hash
  end

  # Allows multiple amplicons to be considered
  # Organises the configuration data by sequence start
  def self.construct_products config_hash
    new_hash = {}
    config_hash[:products].each do |product_name, config|
      new_hash[config[:start]] = config
      new_hash[config[:start]][:output] ||= product_name + ".sam"
      new_hash[config[:start]][:outfile] = File.open(new_hash[config[:start]][:output],'w')
    end
    lengths = new_hash.keys.map!(&:length).uniq
    new_hash[:start_length] = lengths[0]
    warn "Start sequences given must be same length" if lengths.size > 1
    new_hash
  end

  #

end
