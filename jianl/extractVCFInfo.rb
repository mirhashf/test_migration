#ARGV[0] - vcf file
#ARGV[1] - tag
#ARGV[2] - "i"/"f"/"s" (integer/float point/string)
#output to stderr [min, 1st, 2nd, 3rd quarters, max]
if ARGV.size<3
 abort("ERR: arg0-inputs arg1-feature arg2-type")
end

def median(arr)
  lowest = arr.min
  highest = arr.max
  total = arr.inject(:+)
  len = arr.length
  if len > 0
    average = total.to_f / len 
    sorted = arr.sort
    median = len % 2 == 1 ? sorted[len/2] : (sorted[len/2 - 1] + sorted[len/2]).to_f / 2
  else
    median="NA"
  end
  return median
end


ifid=File.open(ARGV[0],"r")
ofid=File.open(ARGV[0]+"."+ARGV[1],"w")
vv=[]
ifid.each{|line|
  if line !~ /\#/
    ll=line.strip.split("\t")
    vv=[]
    ind=ll[8].split(":").index(ARGV[1])
    ofid.write(ll[0]+"\t"+ll[1]+"\t")
    if ind != nil
      ll[9..-1].each{|ii|
        info=ii.split(":")
#        if info[0] != ".\/." && info[0] != "0\/0"
          mm=info[ind]
          if mm != nil
            if ARGV[2] == "i"
              vv.push(mm.to_i)
            elsif ARGV[2] == "f"
              vv.push(mm.to_f)
            
            elsif ARGV[2] == "s"
              vv.push(mm)
            end
          end
#        end
      }      
      ofid.puts(median(vv))
    else
      ofid.puts("NA")
    end

  end
}

ifid.close
ofid.close
