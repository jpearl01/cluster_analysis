#!/usr/bin/env ruby


class Gene
  attr_accessor :sequence, :name, :product, :contig, :start, :end, :organism, :strain
end

gene_list = []

Dir.foreach('/home/josh/Downloads/annotations/') do |f|
  next if f == '.' or f == '..'
  
end
ex_file = '/home/josh/Downloads/annotations/B502.fasta'

