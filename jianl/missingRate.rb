ifid=File.open(ARGV[0],"r")
ofid=File.open(ARGV[0].gsub(".vcf",".missingRate"),"w")

counts=Hash.new(0)

if ARGV.size<1
 abort("ERR: arg0-inputs (has to be vcf)")
end
ifid.each{|line|
  if line !~ /\#/
    ll=line.strip.split("\t")
    gts=ll[9..-1]
    gts.each{|gt|
      counts[gt.split(":")[0]]+=1
    }
    
    ofid.puts("#{ll[0]}\t#{ll[1]}\t#{counts["\.\/\."]/gts.size.to_f}")
    counts.clear

  elsif line =~ /\#/ && line !~ /\#\#/
    ll=line.strip.split("\t")
    ofid.puts(ll[0..1].join("\t")+"\t"+"missingRate")
  end
}

ifid.close
ofid.close
    
