# encoding: utf-8

require 'rubygems'
require 'bundler'
begin
  Bundler.setup(:default, :development)
rescue Bundler::BundlerError => e
  $stderr.puts e.message
  $stderr.puts "Run `bundle install` to install missing gems"
  exit e.status_code
end
require 'rake'

require 'jeweler'
Jeweler::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://guides.rubygems.org/specification-reference/ for more options
  gem.name = "bio-sam-mutation"
  gem.version = '0.4.1'
  gem.homepage = "http://github.com/stveep/bioruby-sam-mutation"
  gem.license = "MIT"
  gem.summary = %Q{Parsing and mutation calling from SAM, CIGAR and MD:Z.}
  gem.description = %Q{Simple classes for parsing SAM, CIGAR and MD:Z strings, including slices. Methods for calling mutations in HGVS format and looking up consequences using Ensembl VEP REST API. Developed for calling mutations at an expected position in an alignment - e.g. Amplicon sequencing of CRISPR-induced mutations.}
  gem.email = "spettitt@gmail.com"
  gem.authors = ["Stephen Pettitt"]
  # dependencies defined in Gemfile
  # gem.required_ruby_version = '>= 1.9.3'
  gem.executables << "mutations"
end
Jeweler::RubygemsDotOrgTasks.new

require 'rake/testtask'
Rake::TestTask.new(:test) do |test|
  test.libs << 'lib' << 'test'
  test.pattern = 'test/**/test_*.rb'
  test.verbose = true
end

desc "Code coverage detail"
task :simplecov do
  ENV['COVERAGE'] = "true"
  Rake::Task['test'].execute
end

task :default => :test

require 'rdoc/task'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "bio-sam-mutation #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
