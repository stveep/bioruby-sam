class Bio::MutationArray < Array
  include VepHgvs
  def to_hgvs(reference_type=nil)
    if length > 0 && reference_type.nil?
      reference_type = first.seqname.match(/^ENS/) ? "c" : "g"
    end
    if length == 1
      first.to_hgvs(reference_type)
    elsif length > 1
      [first.seqname,":",reference_type,".["].join + map(&:to_hgvs).join(";") + "]"
    elsif length == 0
      nil
    end
  end
end
