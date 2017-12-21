# Please require your code below, respecting the naming conventions in the
# bioruby directory tree.
#
# For example, say you have a plugin named bio-plugin, the only uncommented
# line in this file would be
#
#   require 'bio/bio-plugin/plugin'
#
# In this file only require other files. Avoid other source code.

require 'bio'
require 'bio-ensembl-rest'
require 'bio-samtools'
require 'oj'
require 'yaml'
require 'bio-sam-mutation/bio/db/alignment'
require 'bio-sam-mutation/bio/alignment/iterate_pairs'
require 'bio-sam-mutation/bio/alignment/cigar'
require 'bio-sam-mutation/bio/db/tag'
require 'bio-sam-mutation/bio/db/tag/md'
require 'bio-sam-mutation/bio/vephgvs'
require 'bio-sam-mutation/bio/mutation'
require 'bio-sam-mutation/bio/mutation_array'
require 'bio-sam-mutation/bio/mutantallele'

require 'bio-sam-mutation/mutationscli'
