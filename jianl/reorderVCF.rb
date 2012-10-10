if ARGV.size<2 || (ARGV[1] !~ /\.fai/)
 abort("ERR: arg0-inputs arg1-fai file\nUsage:ruby reorderVCF.rb test.vcf test.fai")
end

faifid=File.open(ARGV[1],"r")#("/mnt/scratch0/public/genome/human/hg19.major/hg19.major.fa.fai","r")
ifid=File.open(ARGV[0],"r")
ofid=File.open(ARGV[0].split(".vcf")[0]+".reordered.vcf","w")

vv=Hash.new
ifid.each{|line|
  if line =~ /\#/
    ofid.write(line)
  else
    ll=line.strip.split("\t")
    if vv[ll[0]] == nil
      vv[ll[0]] = []
    end
    vv[ll[0]].push(line)
  end
}

faifid.each{|line|
  ll=line.strip.split("\t")
  if vv[ll[0]] != nil
    vv[ll[0]].each{|rec|
      ofid.write(rec)
    }
  end
}

ifid.close
ofid.close
