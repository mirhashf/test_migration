#RUN AS:
#ruby sepChrFromFasta.rb  <input fasta file> <output path>
##output will be separate contigs in individual files under the output path designated
 
if ARGV.size<2
 abort("ERR: arg0-input fasta file arg1-output path")
end

ifid=File.open(ARGV[0],"r")
outputPath=ARGV[1]
cc=0
ofid=nil
ofn=""

while line=ifid.gets
  if line =~ /\>/
    if cc != 0
      $stderr.puts("#{ofn} written")
      ofid.close
    end
    ofn=outputPath+"\/"+line.strip.gsub("\>","").split(" ")[0]+".fa"
    ofid=File.open(ofn,"w")
    cc+=1
  end
  ofid.write(line)
end

$stderr.puts("#{ofn} written")
ofid.close

ifid.close
