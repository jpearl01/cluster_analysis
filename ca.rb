#!/usr/bin/env ruby

require 'bio'

#################################
#This program is used to take the 'clusterGenes' output from justin's clustering
#perl program and designed to examine the sets of distributed genes.  One of the
#main impetus' for this is to reinclude gene annotation information, like the
#location on the contig the gene is from, and the annotation information 
#provided from RAST (or wherever).
#
#Requirements are the clusterGenes file and the annotation (gff or gbk)
#files for the annotations of each genome.
#################################

#Read in the clusterGenes file

is_new_cluster   = 0
is_gene_seqs     = 0
is_presence_list = 0

File.opt(ARGV[0], "r") do |f|
  while (line = f.gets)
    if line.match(<cluster start>\s+(\S+))
      
    end
  
  end    
end


class Cluster

  attr_accessor :gene_list, :presence_list


end
