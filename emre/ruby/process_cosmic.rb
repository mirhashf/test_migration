#!/usr/bin/ruby
## This script skips lines where the alt alleles are longer than 50 chars
cosmicFile=ARGV[0]

expCols = 10
tabsCount = 0
File.open(cosmicFile, "r").each_line do |line|
  if line.start_with?("#")
    next
  end
  lineTokens = line.split("\t")
  if lineTokens[3].length > 50 || lineTokens[4].length > 5
    next
  end 
  outLine = lineTokens[0..4].join("\t") + "\t"
  tabsCount = 5
  info = lineTokens[7]
  infoTokens = info.split(";")
  i = 0
  infoTokens.each do |it|
    if i > 0
      outLine = outLine + "\t"
      tabsCount += 1
    end
    kv = it.split("=")
    if kv.length == 2
      outLine = outLine + kv[1]
    end
    i += 1
  end
  
  if tabsCount < expCols
    diff = expCols - tabsCount
    outLine.chomp!("\n")
    for j in 1..(diff-1)
      outLine << "\t"
    end
  end
  puts outLine
end
