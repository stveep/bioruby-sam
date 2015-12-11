source "http://rubygems.org"

  gem "bio", ">= 1.4.2"
  # JSON serialisation:
  gem "oj", "~>2.14.0"
  # At the time of writing, the released version 0.2.0 does not include the variation#vep_hgvs method
  # so use this specific commit:
  gem "bio-ensembl-rest", "0.2.0", git: "https://github.com/ALTree/bio-ensembl-rest.git", ref: "c934fa0"

group :development do
  gem "shoulda", ">= 0"
  gem "rdoc", "~> 3.12"
  gem "simplecov", ">= 0"
  gem "jeweler", "~> 2.0.1"
  gem "bundler", "~> 1.10.6"
  gem "test-unit", "~> 3.0.8"
end
