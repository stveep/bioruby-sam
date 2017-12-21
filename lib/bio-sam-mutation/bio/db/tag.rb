Bio::DB::Tag.class_eval do
  def initialize(tag_string=nil)
      set(tag_string) if tag_string
  end
end
