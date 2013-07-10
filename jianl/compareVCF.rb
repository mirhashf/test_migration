#!/usr/bin/env ruby                                                                                            
$VERBOSE = nil

# ##############################################################################                               
# PURPOSE                                                                                                      
# ##############################################################################                               
# Compare multiple vcfs and generate venn diagrams by calling and extracting relavant metrics from GATK VariantEval tool.

# ##############################################################################                               
# REQUIRED LIBRARIES                                                                                           
# ##############################################################################                               
require 'getoptlong'

# ##############################################################################                               
# CONSTANTS                                                                                                    
# ##############################################################################                               

# ##############################################################################                               
# HELPER FUNCTIONS                                                                                             
# ##############################################################################                               
# Process command line args                                                         
def processArguments()
  optsArray = [
                ['--refFasta', '-r', GetoptLong::REQUIRED_ARGUMENT],
                ['--gatkPath', '-g', GetoptLong::REQUIRED_ARGUMENT],
                ['--dbSNP', '-d', GetoptLong::REQUIRED_ARGUMENT],
                ['--vcfFile', '-i', GetoptLong::REQUIRED_ARGUMENT],
                ['--vcfName', '-n', GetoptLong::REQUIRED_ARGUMENT],
                ['--variantClass', '-c', GetoptLong::OPTIONAL_ARGUMENT],
                ['--leftAlign', '-l', GetoptLong::OPTIONAL_ARGUMENT],
                ['--outputDir', '-o', GetoptLong::OPTIONAL_ARGUMENT],
                ['--plotVenn', '-p', GetoptLong::OPTIONAL_ARGUMENT],
                ['--sqrt', '-s', GetoptLong::OPTIONAL_ARGUMENT],
                ['--forceGenFile', '-f', GetoptLong::OPTIONAL_ARGUMENT],
                ['--trueVariantCounts', '-t', GetoptLong::OPTIONAL_ARGUMENT],
                ['--lenient','-u',GetoptLong::OPTIONAL_ARGUMENT],
                ['--help', '-h', GetoptLong::NO_ARGUMENT]
              ]
  progOpts = GetoptLong.new(*optsArray)
  optsHash = Hash.new
  progOpts.each {|opt,arg|
    optsHash[opt]=arg
  }

  usage() if(optsHash.empty? or optsHash.key?('--help') or !optsHash.key?('--refFasta') or !optsHash.key?('--gatkPath') or !optsHash.key?('--dbSNP') or !optsHash.key?('--vcfFile') or !optsHash.key?('--vcfName'))
  return optsHash
end

def mem_usage()
  `ps -p#{$$} -orss`.split[1].to_i
end

def usage(msg='')
  puts "\n#{msg}\n" unless(msg.empty?)
  puts "                                                                                                         
                                                                                                                 
PROGRAM DESCRIPTION:
  Compare multiple vcfs and generate venn diagrams by calling and extracting relavant metrics from GATK VariantEval tool.

  COMMAND LINE ARGUMENTS:                                                                                        
    --refFasta       | -r    => referece fasta file                                                  
    --gatkPath       | -g    => absolute path to GATK's GenomeAnalysisTK.jar
    --dbSNP          | -d    => dbSNP vcf file
    --vcfFile        | -i    => input vcf files for comparison, separated by comma
    --vcfName        | -n    => names for the input vcf files, separated by comma and should be ordered according to the files
    --variantClass   | -c    => specify whether the comparison should be performed for SNPs(\"snp\") or indels(\"indel\") (optional, default is snp)
    --leftAlign      | -l    => if comparing indels, force to apply left alignment (default is true)
    --outputDir      | -o    => output directory (optional, default is the current dir)
    --plotVenn       | -p    => plot the Venn daigrams or not (optional, default is false) 
    --sqrt           | -s    => apply sqrt on the weight to generate the Venn diagrams (only use it if the plotting fails, defaul is false)
    --trueVariantCounts           |-t     => number of true snp or indels in dbSNP (should correspond to --variantClass input, optional)
    --lenient        | -u    => apply LENIENT_VCF_PROCESSING for checking on vcf
    --forceGenFile   | -f    => force to re-generate the snp or indel vcf file for the input, even if the file is existing.
  USAGE:                                                                                                         
  ruby compareVCF.rb -r ref.fa -g pathTogatk/dist -d dbSNP_132.hg19.vcf -i input1.vcf,input2.vcf -n I1,I2 -c snp -o outputPath/ 

  EXAMPLE:

  ruby compareVCF.rb -r /mnt/scratch0/public/genome/human/hg19.major/hg19.major.fa -g /home/jianl/work/seqalto/third-party/gatk/dist -d /mnt/scratch0/public/genome/human/hg19/dbsnp/dbsnp_132.vcf -i /mnt/scratch1/jianl/stanford/cuiping/bina/genotyped.reordered.vcf,/mnt/scratch1/jianl/stanford/cuiping/gatk/TEST-blood-CEU-gatk.reordered.vcf -n BINA,HugeSeq -c snp -o /home/jianl/work/test/

"

  exit(134)
end


def parseGATKtable(ifid,tableName,colName,rods,rodProperty)
  while line=ifid.gets
    if line =~ /\#\:GATK/ && line =~ /#{tableName}/
      break
    end
  end

  tableHeader=ifid.gets.strip.split

  compRodInd=tableHeader.index("CompRod")
  evalRodInd=tableHeader.index("JexlExpression") #("EvalRodInd")                                                                                                            
  noveltyInd=tableHeader.index("Novelty")

  colInd=[]
  colName.each{|col|
    colInd.push(tableHeader.index(col))
    rodProperty.push(col)
  }

  flag=0
  while line=ifid.gets
    ll=line.strip.split
    if ll.size == tableHeader.size && ll[evalRodInd] != "none" && ll[compRodInd] =="dbsnp"
      colInd.each{|ind|
        rods[[ll[compRodInd],ll[evalRodInd]]][ll[noveltyInd]].push(ll[ind])
      }
      if flag == 0
        flag=1
      end
    else
      if ll[evalRodInd] == "none"
        next
      end
      if flag == 1
        break
      end
    end
  end
end

# ------------------------------------------------------------------------------                                 

# ##############################################################################                                 
# MAIN                                                                                                           
# ##############################################################################                                 

$stderr.puts "#{Time.now} BEGIN (Mem: #{mem_usage()})"
optsHash = processArguments()


scale=1#1000000
digitNo=1#10

evalRods=optsHash['--vcfName'].strip.split(",")
evalFiles=optsHash['--vcfFile'].strip.split(",")
gatkPath=optsHash['--gatkPath'].strip
refFasta=optsHash['--refFasta'].strip
dbSNPvcf=optsHash['--dbSNP'].strip
outputPath=(!optsHash.key?('--outputDir')) ? "." : optsHash['--outputDir'].strip
variantClass=(!optsHash.key?('--variantClass')) ? "snp" : optsHash['--variantClass'].strip
plotVennFlag=(!optsHash.key?('--plotVenn')) ? FALSE : optsHash['--plotVenn'].strip
sqrtFlag=(!optsHash.key?('--sqrt')) ? FALSE : optsHash['--sqrt'].strip
genFileFlag=(!optsHash.key?('--forceGenFile')) ? FALSE : optsHash['--forceGenFile'].strip
lenientFlag=(!optsHash.key?('--lenient')) ? FALSE : optsHash['--lenient'].strip
leftAlignFlag=(!optsHash.key?('--leftAlign')) ? TRUE : optsHash['--lenient'].strip

evalCombs=[]
for ii in 1...evalRods.size
  evalRods.combination(ii).to_a.each{|cc|
    evalCombs.push(cc)
  }
end

evalCombs.push(["Intersection"])
selectString=((((["-select 'set==\""]*evalCombs.size).zip(evalCombs.map{|ee| ee.join("-")})).map{|p1| p1.join}.zip((["\"' -selectName "]*evalCombs.size)).map{|p2| p2.join}).zip(evalCombs.map{|ee| ee.join("-")})).map{|p3| p3.join}

novelty=["novel","known"]

rods=Hash.new
compRods=[]
compRodSize=Hash.new

rodProperty=[]

evalSelectedFiles=[]
gatkCMD_u=""
if (lenientFlag)
  gatkCMD_u=" -U LENIENT_VCF_PROCESSING"
end

if variantClass.downcase == "snp"
  featureVenn="nSNPs"
  if !optsHash.key?('--trueVariantCounts')
    compRodSize["dbsnp"]= (`grep -c VC=SN #{dbSNPvcf}`).strip.to_i #28833350 
    $stderr.puts("number of #{variantClass} in #{dbSNPvcf}:\t#{compRodSize["dbsnp"]}")
  else
    compRodSize["dbsnp"]=optsHash['--trueVariantCounts'].strip.to_i
  end
  ii=0
  evalFiles.each{|ff|
    snpFile="#{outputPath}/#{evalRods[ii]}.#{ff.split("\/")[-1]}.snp"
    if File.exists?(snpFile) && !genFileFlag
      $stderr.puts("#{snpFile} exists ... not regenerating")
    else
      system("java -jar #{gatkPath}/GenomeAnalysisTK.jar -T SelectVariants -R #{refFasta} --variant #{ff} -o #{snpFile} -selectType SNP#{gatkCMD_u}")
    end
      evalSelectedFiles.push(snpFile)
    ii+=1
  }

elsif variantClass.downcase == "indel"
  featureVenn="nIndels"
  if !optsHash.key?('--trueVariantCounts')
    compRodSize["dbsnp"]= [(`grep -c VC=INDEL #{dbSNPvcf}`).strip.to_i,(`grep -c VC=DIV #{dbSNPvcf}`).strip.to_i].max
    $stderr.puts("number of #{variantClass} in #{dbSNPvcf}:\t#{compRodSize["dbsnp"]}")
  else
    compRodSize["dbsnp"]=optsHash['--trueVariantCounts'].strip.to_i
  end
  
  ii=0
  evalFiles.each{|ff|
    indelFile="#{outputPath}/#{evalRods[ii]}.#{ff.split("\/")[-1]}.indel"
    indelFile_left="#{outputPath}/#{evalRods[ii]}.#{ff.split("\/")[-1]}.indel.leftAligned"
    if File.exists?(indelFile) && !genFileFlag
      $stderr.puts("#{indelFile} exists ... not regenerating")
    else
      system("java -jar #{gatkPath}/GenomeAnalysisTK.jar -T SelectVariants -R #{refFasta} --variant #{ff} -o #{indelFile} -selectType INDEL#{gatkCMD_u}")
    end
    if(leftAlignFlag)
      system("java -jar #{gatkPath}/GenomeAnalysisTK.jar -T LeftAlignVariants -R #{refFasta} --variant #{indelFile} -o #{indelFile_left}#{gatkCMD_u}")
      indelFile=indelFile_left
    end
  
    evalSelectedFiles.push(indelFile)
    ii+=1
  }
#  evalFiles.map! {|ff| "#{outputPath}/#{ff.split("\/")[-1]}.indel"}
end

combRodFile=(["-V:"]*evalRods.size).zip(evalRods.zip(evalSelectedFiles).map {|rodFile| rodFile.join(" ")}).map {|evalString| evalString.join("")}

combRodFile2=(["--eval:"]*evalRods.size).zip(evalRods.zip(evalSelectedFiles).map {|rodFile| rodFile.join(" ")}).map {|evalString| evalString.join("")}

system("java -jar #{gatkPath}/GenomeAnalysisTK.jar -R #{refFasta} -T CombineVariants  #{combRodFile.join(" ")} -priority #{evalRods.join(",")} -o #{outputPath}/combined.#{evalRods.join("_")}.#{variantClass}.vcf -setKey set#{gatkCMD_u}")
system("java -jar #{gatkPath}/GenomeAnalysisTK.jar -T VariantEval -R #{refFasta} -D #{dbSNPvcf} #{selectString.join(" ")} -o #{outputPath}/combeval.#{evalRods.join("_")}.#{variantClass}.report -eval #{outputPath}/combined.#{evalRods.join("_")}.#{variantClass}.vcf --evalModule GenotypeConcordance -l INFO")

##system("java -jar #{gatkPath}/GenomeAnalysisTK.jar  -T VariantEval -R #{refFasta} -D #{dbSNPvcf} #{combRodFile2.join(" ")} -o weval.#{evalRods.join("_")}.#{variantClass}.report -l INFO")

ifid=File.open("#{outputPath}/combeval.#{evalRods.join("_")}.#{variantClass}.report","r")
statsFid=File.open("#{outputPath}/stats.#{evalRods.join("_")}.#{variantClass}.txt","w")

while line=ifid.gets
  if line =~ /\#\:GATK/ && line =~ /CompOverlap/
    break
  end
end

compOverlapHeader=ifid.gets.strip.split
compRodInd=compOverlapHeader.index("CompRod")  
evalRodInd=compOverlapHeader.index("JexlExpression") #("EvalRodInd")
noveltyInd=compOverlapHeader.index("Novelty")

flag=0
while line=ifid.gets
  ll=line.strip.split
  if ll.size == compOverlapHeader.size && ll[evalRodInd] != "none" && ll[compRodInd] =="dbsnp"
    if rods[[ll[compRodInd],ll[evalRodInd]]] == nil      
      rods[[ll[compRodInd],ll[evalRodInd]]]=Hash.new
    end
    rods[[ll[compRodInd],ll[evalRodInd]]][ll[noveltyInd]]=[]
    compRods.push(ll[compRodInd])

    if flag == 0
      flag=1
    end
  else
    if ll[evalRodInd] == "none"
      next
    end
    if flag == 1
      break
    end
  end
end

compRods.uniq!

parseGATKtable(ifid,"CountVariants",["nSNPs","nInsertions","nDeletions","hetHomRatio"],rods,rodProperty)
rodProperty.push("nIndels")
rods.each{|kk,vv|
  vv.each{|nn,aa|
    aa.push(aa[1].to_i+aa[2].to_i)
  }
}

parseGATKtable(ifid,"TiTvVariantEvaluator",["tiTvRatio"],rods,rodProperty)
ifid.close


if plotVennFlag
  require 'rsruby'
  rInstance=RSRuby.instance
  rInstance.eval_R("suppressMessages(library(\"Vennerable\"))")

  vennCount=Hash.new
  combArr=([0,1]*(evalRods.size+1)).combination(evalRods.size+1).to_a.uniq.sort
  
  weights=Hash.new
  
  compRods.each{|comp|
    if weights[comp] == nil
      weights[comp] = []
    end
  
    knownVector=[]
    
    
    weightByNoveltyString="Weight=c("
    weightAllString="Weight=c("

    for ii in 1...combArr.size
      if ii != combArr.size-1
        if combArr.size > 4 && ii == combArr.size/2-1
          weights[comp].push(rods[[comp,"Intersection"]][novelty[combArr[ii][0]]][rodProperty.index(featureVenn)].to_i) 
        elsif ii == combArr.size/2
          weights[comp].push(compRodSize[comp])
        else
          evalComb=evalRods.values_at(* combArr[ii][1..-1].each_index.select {|jj| combArr[ii][1..-1][jj] ==1})
          weights[comp].push(rods[[comp,evalComb.join("-")]][novelty[combArr[ii][0]]][rodProperty.index(featureVenn)].to_i)
          if ii > combArr.size/2
            knownVector.push(rods[[comp,evalComb.join("-")]]["known"][rodProperty.index(featureVenn)].to_i)
          end
        end
      else
        weights[comp].push(rods[[comp,"Intersection"]][novelty[combArr[ii][0]]][rodProperty.index(featureVenn)].to_i)
        knownVector.push(rods[[comp,"Intersection"]][novelty[combArr[ii][0]]][rodProperty.index(featureVenn)].to_i)
      end
      
    end
    
    weights[comp][combArr.size/2-1] -= rInstance.sum(knownVector)

    weightByNovelty=[]
    weightAll=[]
    
    weights[comp].each{|w| weightByNovelty.push((w.to_f/scale).round(digitNo))}

    for ii in 1...combArr.size
      if sqrtFlag
        weightByNoveltyString+="\""+combArr[ii].join("")+"\" = "+"#{(Math.sqrt(weightByNovelty[ii-1])/scale).round(digitNo)}"
      else
        weightByNoveltyString+="\""+combArr[ii].join("")+"\" = "+"#{weightByNovelty[ii-1]}"
      end
      if ii < combArr.size/2
        weightAll.push(weightByNovelty[ii-1]+weightByNovelty[ii-1+combArr.size/2])
        if sqrtFlag
          weightAllString+="\""+combArr[ii][1..-1].join("")+"\" = "+"#{(Math.sqrt(weightByNovelty[ii-1]+weightByNovelty[ii-1+combArr.size/2])/scale).round(digitNo)}"
        else
          weightAllString+="\""+combArr[ii][1..-1].join("")+"\" = "+"#{weightByNovelty[ii-1]+weightByNovelty[ii-1+combArr.size/2]}"
        end

        if ii != combArr.size/2-1
          weightAllString+=","
        else
          weightAllString+= "\)"
        end
      end
      
      if ii != combArr.size-1
        weightByNoveltyString+=","
      else
        weightByNoveltyString += "\)"
      end
    end  
    rInstance.pdf("#{outputPath}/#{comp}_#{evalRods.join("_")}.#{variantClass}.pdf")
    
    if combArr.size/2 < 10
      $stderr.puts("weightAll:#{weightAll.join(",")}")
      rInstance.eval_R("plot(Venn(SetNames = c(\"#{evalRods.join("\",\"")}\"),#{weightAllString}),doWeights=TRUE,type=\"circles\",show=list(FaceText=\"\",SetLabels=TRUE))") #FaceText=\"weight\"
    end
    
    
    if combArr.size < 10
      rInstance.eval_R("plot(Venn(SetNames = c(\"#{comp}\",\"#{evalRods.join("\",\"")}\"),#{weightByNoveltyString}),doWeights=TRUE,type=\"circles\",show=list(FaceText=\"\",SetLabels=TRUE))")
    else
      rInstance.eval_R("plot(Venn(SetNames = c(\"#{comp}\",\"#{evalRods.join("\",\"")}\"),#{weightByNoveltyString}),doWeights=FALSE,type=\"ellipses\")")
    end
    
    rInstance.eval_R("dev.off()")
    
  }
end

ind=0
rodProperty.each{|prop|
  $stderr.puts("Method\t#{prop}")
  statsFid.puts("Method\t#{prop}")
  $stderr.puts("\tall\tknown\tnovel")
  statsFid.puts("\tall\tknown\tnovel")
  rods.each{|kk,vv|
    if kk[0] == "dbsnp"
      $stderr.puts("#{kk[1]}\t#{vv["all"][ind]}\t#{vv["known"][ind]}\t#{vv["novel"][ind]}")
      statsFid.puts("#{kk[1]}\t#{vv["all"][ind]}\t#{vv["known"][ind]}\t#{vv["novel"][ind]}")
    end
  }
  ind+=1
  $stderr.puts("")
  statsFid.puts("")
}
statsFid.close

$stderr.puts "#{Time.now} DONE"
exit(0)
