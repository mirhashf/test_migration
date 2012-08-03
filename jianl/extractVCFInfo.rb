#ARGV[0] - vcf file
#ARGV[1] - tag
#ARGV[2] - "i"/"f" (integer/float point)
#output to stderr [min, 1st, 2nd, 3rd quarters, max]

ifid=File.open(ARGV[0],"r")
vv=[]
ifid.each{|line|
  if line !~ /\#/
    if ARGV[2] == "i"
      mm=/#{ARGV[1]}=(\d+)/.match(line)
      if mm != nil
        vv.push(mm[1].to_i)
      end
    elsif ARGV[2] == "f"
      mm=/#{ARGV[1]}=([-+]?[0-9]*\.?[0-9]+)/.match(line)
      if mm != nil
        vv.push(mm[1].to_f)
      end
    end

  end
}

# $stderr.puts(vv.size)
len=vv.size
vv.sort!
$stderr.puts("[#{vv[0]}, #{vv[len/4]}, #{vv[len/2]}, #{vv[len*3/4]}, #{vv[-1]}]")
