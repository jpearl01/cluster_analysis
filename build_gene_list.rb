#!/usr/bin/env ruby


class Gene
  attr_accessor :sequence, :name, :product, :contig, :start, :end, :organism, :strain
end

gene_list = []

Dir.foreach('/home/josh/Downloads/annotations/') do |f|
  next if f == '.' or f == '..'
  
end
ex_file = '/home/josh/Downloads/annotations/B502.fasta'

fh = File.open(ex_file, 'r')

genome_name = File.basename(ex_file, '.*')

fh.each_line do |l|
  
  if l.match('>')
    next if l.match(' Contig ')
    gene = Gene.new
    l.match('>([\S]+)')
    gene.name = $1
    gene.contig = gene.name.split('_')[0]
    gene.name.gsub!(gene.contig, genome_name)
    puts gene.name
    l.match('(product)=[\'\"]{2}(.+)[\'\"]{2}')
    gene.product = $2
    puts gene.product
  end

end

