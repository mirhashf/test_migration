if ARGV.size<2 || (ARGV[0] =~ /\*/ && ARGV[0] !~ /\\*/)
 abort("ERR: arg0-inputs arg1-output")
end
system("ls #{ARGV[0]}>file.ls")
ifid=File.open("file.ls","r")
inputs=""
ifid.each{|line|
  inputs+="I=#{line.strip} "
}
ifid.close
$stderr.puts("java -Xms15g -Xmx15g -jar ~/downloads/picard/dist/MergeSamFiles.jar #{inputs} OUTPUT=#{ARGV[1]} SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT USE_THREADING=T MSD=true")
system("java -Xms15g -Xmx15g -jar ~/downloads/picard/dist/MergeSamFiles.jar #{inputs} OUTPUT=#{ARGV[1]} SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT USE_THREADING=T MSD=true")
