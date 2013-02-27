#!/usr/bin/env ruby


class Gene
  #Here the db_xref is an array
  attr_accessor :sequence, :name, :product, :contig, :start, :end, :organism, :strain, :db_xref, :is_complement
end

gene_list = []

Dir.foreach('data/annotations/') do |f|
  next if f == '.' or f == '..'
  
end
ex_file = 'data/annotations/B502.fasta'

fh = File.open(ex_file, 'r')

genome_name = File.basename(ex_file, '.*')

fh.each_line do |l|
  
  if l.match('>')
    next if l.match(' Contig ')
    gene = Gene.new
    gene.strain = genome_name
    l.match('>([\S]+)')
    gene.name = $1
    gene.contig = gene.name.split('_')[0]
    gene.name.gsub!(gene.contig, genome_name)
#    puts gene.name
    l.match('(product)=[\'\"]{2}(.+)[\'\"]{2}')
    gene.product = $2
#    puts gene.product
    if l.match('complement')
      gene.is_complement = true
    else
      gene.is_complement = false
    end
    l.match('loc=[complement]*\(*(\d+)\.\.(\d+)\)*')
    gene.start = $1.to_i
    gene.end   = $2.to_i
#    puts "Gene start " + gene.start.to_s + " gene end "+ gene.end.to_s
    #Adding a check to see if more than one db_xref
    gene.db_xref = l.scan(/db_xref="[^"]+/)
    gene.db_xref.each_entry {|s| s.gsub!('db_xref="', '')}
#    puts "The list of db_xrefs is " + gene.db_xref.join(' ')
    gene_list.push(gene)
  end

end

puts "Total genes are: "+ gene_list.size.to_s
puts gene_list.first.inspect
