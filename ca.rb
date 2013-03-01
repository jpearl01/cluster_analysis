#!/usr/bin/env ruby

require 'axlsx'
require 'trollop'

#################################
#This program is used to take the 'clusterGenes' output from justin's clustering
#perl program and designed to examine the sets of distributed genes.  One of the
#main impetus' for this is to reinclude gene annotation information, like the
#location on the contig the gene is from, and the annotation information 
#provided from RAST (or wherever).
#
#Requirements are the clusterGenes file and the directory where annotation 
#files for the annotations of each genome.
#################################

d = '/home/josh/workspace/bioruby/cluster_analysis/data/annotations/'

#Setup options hash
opts = Trollop::options do
  opt :clust_genes, "The mandatory clusterGenes file", type: :string, short: '-c'                     
  opt :directory, "Mandatory directory where the fasta annotation files are (directory cannot have other files)", type: :string, short: '-d'
end

abort("Must have a clusterGenes file defined '-h' or '--help' for usage") if opts[:clust_genes].nil?
abort("The clusterGenes file must actually exist '-h' or '--help' for usage") unless File.exist?(opts[:clust_genes])
abort("Must have a directory with annotations defined '-h' or '--help' for usage") if opts[:directory].nil?
abort("The directory must actually exist '-h' or '--help' for usage") unless Dir.exist?(opts[:directory])

#Define the Cluster and Gene classes

class Cluster
  attr_accessor :gene_list, :presence_list, :name
end

class Gene
  #Here the db_xref is an array
  attr_accessor :sequence, :name, :product, :contig, :start, :end, :organism, :strain, :db_xref, :is_complement
end

is_gene_seqs     = false
is_presence_list = false
all_clusters     = []
gene_hash        = {}


#Read in the clusterGenes file

File.open(opts[:clust_genes], "r") do |f|
  while (line = f.gets)
    if line.match('<cluster start>\s+(\S+)')

      #Initialize the cluster object
      clust               = Cluster.new
      clust.name          = $1
      clust.gene_list     = []
      clust.presence_list = {}

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



#Assuming that in the directory there are only annotation files, and only the number of annotation files for this project
#also have to account for '.' and '..' so minus 2
NUM_STRAINS = Dir.entries(opts[:directory]).size - 2
puts "The number of strains are: #{NUM_STRAINS}"

Dir.foreach(opts[:directory]) do |f|
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
      #This ruby's version of a try/catch block
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


p = Axlsx::Package.new
wb = p.workbook

#Build the sheets needed for this project into the workbook
(1..NUM_STRAINS).each_entry do |i|
  wb.add_worksheet(:name => "cluster_size_#{i}") do |sheet|
    sheet.add_row ["All clusters of size #{i}"]
    sheet.add_row %w(Strain Gene_Name Contig Start End Complement Product)
  end
end

#There is a chance we'll run into clusters that have more genes than the number of strains (duplicate genes and whatnot)
#Adding a worksheet for that eventuality:
wb.add_worksheet(name: "cluster_size_gt_#{NUM_STRAINS}" ) do |sheet|
  sheet.add_row ["All clusters greater than size #{NUM_STRAINS}"]
  sheet.add_row %w(Strain Gene_Name Contig Start End Complement Product)
end

background = ''
#Try adding styles here?
wb.styles do |s|
  background = s.add_style bg_color: 'FFFFCC'
end

#Now populate the worksheets with the different clusters
all_clusters.each_entry do |c|
  c_size = c.gene_list.size
  if c_size <= NUM_STRAINS
    wb.sheet_by_name("cluster_size_#{c_size}").add_row [c.name]
    c.gene_list.each_entry do |g|
      wb.sheet_by_name("cluster_size_#{c_size}").add_row [ gene_hash[g].strain, 
                                                           gene_hash[g].name, 
                                                           gene_hash[g].contig, 
                                                           gene_hash[g].start, 
                                                           gene_hash[g].end, 
                                                           gene_hash[g].is_complement, 
                                                           gene_hash[g].product ], style: background
    end
  else
    wb.sheet_by_name("cluster_size_gt_#{NUM_STRAINS}").add_row [c.name, "Total Size: #{c_size}"]
    c.gene_list.each_entry do |g|
      wb.sheet_by_name("cluster_size_gt_#{NUM_STRAINS}").add_row [ gene_hash[g].strain, 
                                                                   gene_hash[g].name, 
                                                                   gene_hash[g].contig, 
                                                                   gene_hash[g].start, 
                                                                   gene_hash[g].end, 
                                                                   gene_hash[g].is_complement, 
                                                                   gene_hash[g].product ]
    end
  end
end

p.serialize("trial_gene_list.xlsx")
