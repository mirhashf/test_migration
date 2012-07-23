#RUN AS:
#ruby sepChrFromFasta.rb  <input fasta file>
##output will be separate contigs in individual files under the same directory as input
 
ifid=File.open(ARGV[0],"r")
cc=0
ofid=nil
ofn=""

while line=ifid.gets
  if line =~ /\>/
    if cc != 0
      $stderr.puts("#{ofn} written")
      ofid.close
    end
    ofn=ARGV[0].split("\/")[0..-2].join("\/")+"\/"+line.strip.gsub("\>","")+".fa"
    ofid=File.open(ofn,"w")
    cc+=1
  end
  ofid.write(line)
end

$stderr.puts("#{ofn} written")
ofid.close

ifid.close
