ifid=File.open(ARGV[0],"r")
ofid=File.open(ARGV[0].split("\/")[-1].gsub("hg19","human_g1k_v37_decoy"),"w")
col=ARGV[1].to_i

corrfid=File.open("/mnt/scratch0/public/genome/human//map_hg19_decoy.txt","r")
corrfid.gets
corr=Hash.new
corrfid.each{|line|
  ll=line.strip.split("\t")
  corr[ll[1]] = ll[2]
}
corrfid.close

ifid.each{|line|
  if ARGV[1] != nil
    ll=line.strip.split("\t")
    ll[col]=corr[ll[col]]  
    ofid.puts(ll.join("\t"))
  else
    chr=/chr(\w+)/.match(line)
    if chr[0] != nil
      ofid.write(line.gsub(chr[0],corr[chr[0]]))
    else
      $stderr.write(line)
    end
  end
    
}

ifid.close
ofid.close
