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

#Define the Cluster and Gene classes

class Cluster
  attr_accessor :gene_list, :presence_list, :name
end


class Gene
  attr_accessor :sequence, :name, :product, :contig, :start, :end, :organism, :strain
end

#Read in the clusterGenes file

is_gene_seqs     = false
is_presence_list = false
all_clusters     = Array.new

File.open(ARGV[0], "r") do |f|
  while (line = f.gets)
    if line.match('<cluster start>\s+(\S+)')
      #Initialize the cluster object
      clust = Cluster.new
      clust.name=$1
      clust.gene_list = Array.new
      clust.presence_list = Hash.new
      all_clusters.push(clust)
    end

    #We are in the part of the clusterGenes file that lists the different included genes fasta style
    is_gene_seq = true if line.match('<genes start>')
    is_gene_seq = false if line.match('<genes stop>')

  end  
end

puts all_clusters.size
