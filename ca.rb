#!/usr/bin/env ruby

require 'axlsx'
require 'trollop'
require 'yaml'

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
  opt :clust_genes, "Mandatory clusterGenes file (output from Justin's clsuterGenes.pl program)", type: :string, short: '-c'                     
  opt :directory, "Mandatory directory where the fasta annotation files are (directory cannot have other files)", type: :string, short: '-d'
  opt :gene_map, "Optional gene map file, if you need to replace gene names.  Of the form:\ncurrentName_1\tnewName_1\ncurrentName_2\tnewName_2\n", type: :string, short: '-m'
  opt :genome_groups, "Optional clade grouping file in YAML format.", type: :string, short: '-g'
  opt :custom_list, "Optional list of genes, these clusters of interest will be put in the \"Custom_list\" worksheet.", type: :string, short: '-l'
end

abort("Must have a clusterGenes file defined '-h' or '--help' for usage") if opts[:clust_genes].nil?
abort("The clusterGenes file must actually exist '-h' or '--help' for usage") unless File.exist?(opts[:clust_genes])
abort("Must have a directory with annotations defined '-h' or '--help' for usage") if opts[:directory].nil?
abort("The fasta annotation directory must actually exist '-h' or '--help' for usage") unless Dir.exist?(opts[:directory])
if !opts[:custom_list].nil?
  abort("You provided a custom_list but it doesn't exist!") unless File.exists?(opts[:custom_list])
end

#Define the Cluster and Gene classes

class Cluster
  attr_accessor :gene_list, :presence_list, :name
end

class Gene
  #Here the db_xref is an array
  attr_accessor :sequence, :name, :product, :contig, :start, :end, :organism, :strain, :db_xref, :is_complement, :size, :clade
end

#Initialize a hash with a key for each gene class (1..n) where n is the number of strains in the project, and each gene class has an entry for each group, and if it is cross-clade
#Input is the genome_groups, and the list of Genomes
def init_gene_class_groups(g_groups={}, g_list=[])
  c_class = {}
  Array(1..g_list.size+1).each_entry do |e|
    c_class[e] = {}
    g_groups.keys.each_entry do |k|
      c_class[e][k] = 0
    end
    c_class[e]['cross_clade'] = 0
  end
  return c_class
end


#Method to track the number of clusters from either a single clade, or cross clade
def track_groups(c_class_hash, cluster)

end

is_gene_seqs     = false
is_presence_list = false
all_clusters     = []
genome_list      = []
gene_hash        = {}
gene_map         = {}     
genome_groups    = {}

#Read in the custom genes list, if exists
custom_genes = IO.read(opts[:custom_list]).split unless opts[:custom_list].nil?


#Read in the map file
if !opts[:gene_map].nil?
  File.open(opts[:gene_map], "r") do |f|
    while (line = f.gets)
      a = line.split
      gene_map[a[0]] = a[1]
    end
  end
end

#Read in the grouping file
genome_groups  = YAML.load_file(opts[:genome_groups]) if !opts[:genome_groups].nil?
genome_groups.each { |key, arr|
  arr.map! {|v| v.to_s}
}


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
      if gene_map.has_key?($1)
        clust.gene_list.push(gene_map[$1])
      else
        clust.gene_list.push($1)
      end
    end

    #We are in the part of the clusterGenes file that lists the different matched strains
    is_presence_list = true if line.match('<strains start>')
    is_presence_list = false if line.match('<strains stop>')

    #populate the presence hash
    if is_presence_list 
      clust.presence_list[line.split[0]] = line.split[1] unless line.match('<strains start>')
    end
  end
end

puts "Total number of clusters is: " + all_clusters.size.to_s



#Assuming that in the directory there are only annotation files, and only the number of annotation files for this project
#also have to account for '.' and '..' so minus 2
NUM_STRAINS = Dir.entries(opts[:directory]).size - 2
puts "The number of strains are: #{NUM_STRAINS}"

#Read in all of our annotations, put them in the gene_hash
Dir.foreach(opts[:directory]) do |f|
  next if f == '.' or f == '..'

  fh = File.open(d+f, 'r')
  genome_name = File.basename(f, '.*')
  genome_list.push(genome_name)
  puts "#{f} is the current file"

  fh.each_line do |l|
    
    if l.match('>')

      #just incase we missed a contig entry, though that is unlikely.
      next if l.match(' Contig ')
      gene = Gene.new

      #parse out the current gene name, update it, and set which contig the gene came from
      #This is a little tricky... have to pay attention to the gene map names, they should match the genome name
      gene.strain = genome_name
      l.match('>([\S]+)')
      gene.name = $1
      gene.contig = gene.name.split('_')[0]
      gene.name.gsub!(gene.contig, genome_name)
      if gene_map.has_key?(gene.name)
        gene.name = gene_map[gene.name]
      end

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
      $stderr.puts "There is probably a complex gene with name #{gene.name}, currently have to manually fix this" if l.match('join\(')
      gene.size = (gene.end - gene.start).abs

      #Adding a check to see if more than one db_xref
      gene.db_xref = l.scan(/db_xref="[^"]+/)
      gene.db_xref.each_entry {|s| s.gsub!('db_xref="', '')}

      #Finally add this gene to the gene hash
      gene_hash[gene.name]=gene
    end
    
  end

end

#Initialize the gene classes based on the yaml file
init_gene_class_groups(genome_groups, genome_list) unless opts[:genome_groups].nil? 

#Add the first workbook, this holds the original annotation data for each cluster
p = Axlsx::Package.new
wb = p.workbook


#Build the sheets needed for this project into the workbook
(1..NUM_STRAINS).each_entry do |i|

  wb.add_worksheet(:name => "cluster_size_#{i}") do |sheet|
    sheet.add_row ["All clusters of size #{i}"]
    sheet.add_row %w(Strain Gene_Name Contig Start End Size Complement Product)
  end

end

#On second thought, really all we need is a new sheet in the workbook which has the presence/absence lists
wb.add_worksheet(name: "Presence_Absence") do |sheet|
  sheet.add_row ["Cluster Presence/Absence"]
  row_titles = genome_list.sort.unshift('cluster_names')
  sheet.add_row row_titles
end

#Need to add in another worksheet which will keep track of the names of the genes present in the cluster
wb.add_worksheet(name: "Pres_abs_names") do |sheet|
  sheet.add_row ["Cluster Presence/Absence"]
  row_titles = genome_list.sort.unshift('cluster_names')
  sheet.add_row row_titles
end

#There is a chance we'll run into clusters that have more genes than the number of strains (duplicate genes and whatnot)
#Adding a worksheet for that eventuality, this isn't an issue in the presence workbook, as if it is there at least once it is a 1, 0 otherwise
wb.add_worksheet(name: "cluster_size_gt_#{NUM_STRAINS}" ) do |sheet|
  sheet.add_row ["All clusters greater than size #{NUM_STRAINS}"]
  sheet.add_row %w(Strain Gene_Name Contig Start End Size Complement Product)
end

#And add a custom sheet for when someone wants a list of clusters containing specific genes
wb.add_worksheet(name: "Custom_list") do |sheet| 
  sheet.add_row ["Custom cluster list"]
  sheet.add_row %w(Strain Gene_Name Contig Start End Size Complement Product)
end

background = ''
group1     = ''
group2     = ''
#Try adding styles here?
wb.styles do |s|
  background = s.add_style bg_color: 'FFFFCC'
  group1     = s.add_style bg_color: '87CEFA'
  group2     = s.add_style bg_color: 'ecb0b0'
end

#Need a hash to contain gene classes with more than the normal core max (duplicates, broken gene annotations etc.) so I can print out some stats on it
#Also want to find the number of clade specific genes, and the number of genes that are cross clade
gene_class_gt_core = {}
clade_specific     = {}

#Now populate the worksheets with the different clusters
all_clusters.each_entry do |c|
  c_size = c.gene_list.size
  #Annoyingly, I'm not sure how to do the pres/absence thing except by a little bit of a wonky enumerable
  pres = []
  pres.push(c.name)
  c.presence_list.sort_by {|genome, pres| genome}.each_entry do |p|
    pres.push(p[1])
  end
  wb.sheet_by_name("Presence_Absence").add_row pres

  #ok, gotta check each clade's members
  g1 = 0
  g2 = 0

  c.presence_list.each_key do |h|
    if c.presence_list[h] == "1"
      h = h.to_s unless h.class == String
      g2 += 1 if genome_groups['seroresistant'].include? h
      g1 += 1 if genome_groups['serosensitive'].include? h
    end
  end

  if clade_specific[c_size].nil?
    clade_specific[c_size] = {'seroresistant' => 0, 'serosensitive' => 0, 'cross_clade' => 0}
  end

  if g2 == 0 && g1 > 0
    clade_specific[c_size]['serosensitive'] += 1
  elsif g2 > 0 && g1 == 0
    clade_specific[c_size]['seroresistant'] += 1
  else
    clade_specific[c_size]['cross_clade']   += 1
  end
  
  is_custom_gene = false

  if c_size <= NUM_STRAINS

    wb.sheet_by_name("cluster_size_#{c_size}").add_row [c.name]
    c.gene_list.each_entry do |g|
      abort("#{g.name}") if gene_hash[g].nil?
      #TODO: Make this be able to read in arbitrary genome groupings
      name = gene_hash[g].strain.to_s unless gene_hash[g].strain.class == String
      name = gene_hash[g].strain if gene_hash[g].strain.class == String
      bkgrnd = 0
      if genome_groups['seroresistant'].include? name
        bkgrnd = group2
      elsif genome_groups['serosensitive'].include? name
        bkgrnd = group1
      else
        bkgrnd = background
      end

      #check if this is a custom gene
      is_custom_gene = true if !custom_genes.nil? && (custom_genes.include? gene_hash[g].name)

      wb.sheet_by_name("cluster_size_#{c_size}").add_row [ gene_hash[g].strain, 
                                                           gene_hash[g].name, 
                                                           gene_hash[g].contig, 
                                                           gene_hash[g].start, 
                                                           gene_hash[g].end,
                                                           gene_hash[g].size,
                                                           gene_hash[g].is_complement, 
                                                           gene_hash[g].product ], style: bkgrnd
    end
  else
    #initialize the gene class size to zero, if it doesn't exist
    gene_class_gt_core[c_size] = 0 if gene_class_gt_core[c_size].nil?
    gene_class_gt_core[c_size] += 1

    wb.sheet_by_name("cluster_size_gt_#{NUM_STRAINS}").add_row [c.name, "Total Size: #{c_size}"]
    c.gene_list.each_entry do |g|
      
      #TODO: Make this be able to read in arbitrary genome groupings
      name = gene_hash[g].strain.to_s unless gene_hash[g].strain.class == String
      name = gene_hash[g].strain if gene_hash[g].strain.class == String
      bkgrnd = 0
      if genome_groups['seroresistant'].include? name
        bkgrnd = group2
      elsif genome_groups['serosensitive'].include? name
        bkgrnd = group1
      else 
        bkgrnd = background
      end

      #check if this is a custom gene
      is_custom_gene = true if !custom_genes.nil? && (custom_genes.include? gene_hash[g].name)

      wb.sheet_by_name("cluster_size_gt_#{NUM_STRAINS}").add_row [ gene_hash[g].strain, 
                                                                   gene_hash[g].name, 
                                                                   gene_hash[g].contig, 
                                                                   gene_hash[g].start, 
                                                                   gene_hash[g].end,
                                                                   gene_hash[g].size,
                                                                   gene_hash[g].is_complement, 
                                                                   gene_hash[g].product ], style: bkgrnd
    end
  end
  if is_custom_gene
    wb.sheet_by_name("Custom_list").add_row [c.name, "Total Size: #{c_size}"]
    c.gene_list.each_entry do |g|
      #TODO: Make this be able to read in arbitrary genome groupings
      name = gene_hash[g].strain.to_s unless gene_hash[g].strain.class == String
      name = gene_hash[g].strain if gene_hash[g].strain.class == String
      bkgrnd = 0
      if genome_groups['seroresistant'].include? name
        bkgrnd = group2
      elsif genome_groups['serosensitive'].include? name
        bkgrnd = group1
      else 
        bkgrnd = background
      end
      wb.sheet_by_name("Custom_list").add_row [ gene_hash[g].strain, 
                                                                   gene_hash[g].name, 
                                                                   gene_hash[g].contig, 
                                                                   gene_hash[g].start, 
                                                                   gene_hash[g].end,
                                                                   gene_hash[g].size,
                                                                   gene_hash[g].is_complement, 
                                                                   gene_hash[g].product ], style: bkgrnd
      
    end
  end
end
#puts clade_specific.to_yaml
#abort("testing")

#print out stats on gene classes with gt core genes
gene_class_gt_core.each do |key, value|
  puts "#{key}\t#{value}"
end

clade_specific.each_entry do |k, v|
  s = v['serosensitive']
  r = v['seroresistant']
  x = v['cross_clade']
  puts k.to_s+"\t"+x.to_s+"\t"+s.to_s+"\t"+r.to_s
end

puts "total clusters seen:"+ all_clusters.size.to_s
p.serialize("clustered_gene_list.xlsx")

