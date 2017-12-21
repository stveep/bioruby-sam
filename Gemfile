source "http://rubygems.org"

  gem "bio", ">= 1.4.2"
  gem "bio-samtools", "~> 2.4" 
  # JSON serialisation:
  gem "oj", "~>2.14"
  # At the time of writing, the released version 0.2.0 does not include the variation#vep_hgvs method
  # so use this specific commit:
  gem "bio-ensembl-rest", "0.2.0", git: "https://github.com/ALTree/bio-ensembl-rest.git", ref: "c934fa0"
  gem "trollop"
  gem "rake"
  gem "net-ping"
  gem "roo", "~>2.5.0"

group :development do
  gem "shoulda", ">= 0"
  gem "rdoc", "~> 3.12"
  gem "simplecov", ">= 0"
  gem "jeweler", "~> 2.0"
  gem "bundler"
  gem "test-unit", "~> 3.0"
  gem "rspec"
end
