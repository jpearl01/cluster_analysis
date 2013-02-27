#!/usr/bin/env ruby

require 'axlsx'

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
  #Here the db_xref is an array
  attr_accessor :sequence, :name, :product, :contig, :start, :end, :organism, :strain, :db_xref, :is_complement
end

#Read in the clusterGenes file

is_gene_seqs     = false
is_presence_list = false
all_clusters     = Array.new
gene_hash        = Hash.new

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
    is_gene_seq = true  if line.match('<genes start>')
    is_gene_seq = false if line.match('<genes stop>')

    #populate the gene list
    if is_gene_seq && line.match('^>(\S+)')
      clust.gene_list.push($1)
    end

    #We are in the part of the clusterGenes file that lists the different matched strains
    is_presence_list = true if line.match('<strains start>')
    is_presence_list = false if line.match('<strains stop>')

    #populate the presence hash
    if is_presence_list 
      clust.presence_list[line.split[0]] = line.split[1]
    end

  end  
end

puts all_clusters.size

d = '/home/josh/workspace/bioruby/cluster_analysis/data/annotations/'

Dir.foreach(d) do |f|
  next if f == '.' or f == '..'

  fh = File.open(d+f, 'r')
  genome_name = File.basename(f, '.*')
  puts "#{f} is the current file"

  fh.each_line do |l|
    
    if l.match('>')

      #just incase we missed a contig entry, unlikely.
      next if l.match(' Contig ')
      gene = Gene.new

      #parse out the current gene name, update it, and set which contig the gene came from
      gene.strain = genome_name
      l.match('>([\S]+)')
      gene.name = $1
      gene.contig = gene.name.split('_')[0]
      gene.name.gsub!(gene.contig, genome_name)

      #parse out and set the product
      l.match('(product)=[\'\"]{1,3}([^\"]+)[\'\"]{1,3}[^;]*')
      gene.product = $2

      #check and set if the gene is a complement
      if l.match('complement')
        gene.is_complement = true
      else
        gene.is_complement = false
      end

      #get and set start and end locations
      l.match('loc=[complement]*\(*(\d+)\.\.(\d+)\)*')
      gene.start = $1.to_i
      gene.end   = $2.to_i

      #Adding a check to see if more than one db_xref
      gene.db_xref = l.scan(/db_xref="[^"]+/)
      gene.db_xref.each_entry {|s| s.gsub!('db_xref="', '')}

      #Finally add this gene to the gene hash
      gene_hash[gene.name]=gene
    end
    
  end

end

#Time to add in the final loop that will print out all the annotations for each cluster, separated by line with only the cluster name.
all_clusters.each_entry do |c|
  puts c.name
  if c.gene_list.size > 0
    c.gene_list.each_entry do |g|
#      if gene_hash.include? g
      begin
        puts '>'+gene_hash[g].name + "\t" + gene_hash[g].product
      rescue
        $stderr.puts "The gene #{g} from cluster #{c.name} does not exist in the gene hash"
        $stderr.puts gene_hash[g].inspect
        next
      end
    end
  else 
    $stderr.puts "The cluster #{c.name} does not have any genes in it"
  end
end
