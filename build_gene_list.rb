#!/usr/bin/env ruby

require 'axlsx'

class Gene
  #Here the db_xref is an array
  attr_accessor :sequence, :name, :product, :contig, :start, :end, :organism, :strain, :db_xref, :is_complement
end
gene_hash = Hash.new

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
      l.match('(product)=[\'\"]{2}(.+)[\'\"]{2}')
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
break
#  puts gene_hash[genome_name+'_1'].inspect
end

#puts gene_hash['MX1_1'].inspect



p = Axlsx::Package.new

wb = p.workbook
#Build in the sheets needed for this project
(1..31).each_entry do |i|
  wb.add_worksheet(:name => "cluster_size_#{i}") do | sheet |
    sheet.add_row ['some example']
#    gene_hash.each_entry do |g|
#      sheet.add_row [ g[1].strain, g[1].name, g[1].start, g[1].end, g[1].is_complement, g[1].product ]
#    end
  end
end

wb.sheet_by_name("cluster_size_1").add_row ['stuff','does this get added, I hope?']

# do |sheet|
#  puts 'even happen?'
#  sheet
#end
#puts 'wth'

p.serialize("trial_gene_list.xlsx")

