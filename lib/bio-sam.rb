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
require 'bio-sam/sam'
require 'bio-sam/iterate_pairs'
require 'bio-sam/cigar'
require 'bio-sam/mdtag'
require 'bio-sam/mutation'
