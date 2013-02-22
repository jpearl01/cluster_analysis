#!/usr/bin/env ruby

require 'bio'
require 'nokogiri'

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
cg_file = ARGV[0]
doc = Nokogiri::XML(open(cg_file))
