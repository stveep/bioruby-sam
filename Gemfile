source "http://rubygems.org"

  gem "bio", ">= 1.4.2"
  # Using edge version due to a problem with ruby >2.1 in biogems version at time of writing
  gem "bio-samtools", "~>2.3.4", git: "https://github.com/helios/bioruby-samtools.git", ref: "2e77274"
  # JSON serialisation:
  gem "oj", "~>2.14"
  # At the time of writing, the released version 0.2.0 does not include the variation#vep_hgvs method
  # so use this specific commit:
  gem "bio-ensembl-rest", "0.2.0", git: "https://github.com/ALTree/bio-ensembl-rest.git", ref: "c934fa0"
  gem "trollop"
  gem "rake", "~>0.9"

group :development do
  gem "shoulda", ">= 0"
  gem "rdoc", "~> 3.12"
  gem "simplecov", ">= 0"
  gem "jeweler", "~> 2.0"
  gem "bundler", "~> 1.7"
  gem "test-unit", "~> 3.0"
end
