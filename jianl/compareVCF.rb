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
require 'rsruby'
require '/home/jianl/work/script/findoverlap.rd.rb'

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
                ['--outputDir', '-o', GetoptLong::OPTIONAL_ARGUMENT],
                ['--plotVenn', '-p', GetoptLong::OPTIONAL_ARGUMENT],
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
    --outputDir      | -o    => output directory (optional, default is the current dir)
    --plotVenn       | -p    => plot the Venn daigrams or not (optional, default is false) 

  USAGE:                                                                                                         
  ruby compareVCF.rb -r ref.fa -g pathTogatk/dist -d dbSNP_132.hg19.vcf -i input1.vcf,input2.vcf -n I1,I2 -c snp -o outputPath/ 

  EXAMPLE:

  ruby compareVCF.rb -r /mnt/scratch0/public/genome/human/hg19.major/hg19.major.fa -g /home/jianl/work/gatk_jun12/seqalto/third-party/gatk/dist -d /mnt/scratch0/public/genome/human/hg19/dbsnp/dbsnp_132.vcf -i /mnt/scratch1/jianl/stanford/cuiping/bina/genotyped.reordered.vcf,/mnt/scratch1/jianl/stanford/cuiping/gatk/TEST-blood-CEU-gatk.reordered.vcf -n BINA,GATK -c snp -o /home/jianl/work/test/ -p true

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
digitNo=10

evalRods=optsHash['--vcfName'].strip.split(",")
evalFiles=optsHash['--vcfFile'].strip.split(",")
gatkPath=optsHash['--gatkPath'].strip
refFasta=optsHash['--refFasta'].strip
dbSNPvcf=optsHash['--dbSNP'].strip
outputPath=(!optsHash.key?('--outputDir')) ? "." : optsHash['--outputDir'].strip
variantClass=(!optsHash.key?('--variantClass')) ? "snp" : optsHash['--variantClass'].strip
plotVennFlag=(!optsHash.key?('--plotVenn')) ? FALSE : optsHash['--plotVenn'].strip


evalCombs=[]
for ii in 1...evalRods.size
  evalRods.combination(ii).to_a.each{|cc|
    evalCombs.push(cc)
  }
end

evalCombs.push(["Intersection"])
selectString=((((["-select 'set==\""]*evalCombs.size).zip(evalCombs.map{|ee| ee.join("-")})).map{|p1| p1.join}.zip((["\"' -selectName "]*evalCombs.size)).map{|p2| p2.join}).zip(evalCombs.map{|ee| ee.join("-")})).map{|p3| p3.join}

novelty=["novel","known"]
rInstance=RSRuby.instance
rInstance.eval_R("suppressMessages(library(\"Vennerable\"))")

rods=Hash.new
compRods=[]
compRodSize=Hash.new

rodProperty=[]

if variantClass.downcase == "snp"
  featureVenn="nSNPs"
  compRodSize["dbsnp"]= (`grep -c VC=SN #{dbSNPvcf}`).strip.to_i #28833350 
  evalFiles.each{|ff|
    system("java -jar #{gatkPath}/GenomeAnalysisTK.jar -T SelectVariants -R #{refFasta} --variant #{ff} -o #{outputPath}/#{ff.split("\/")[-1]}.snp -selectType SNP")
  }
  evalFiles.map! {|ff| "#{outputPath}/#{ff.split("\/")[-1]}.snp"}
elsif variantClass.downcase == "indel"
  featureVenn="nIndels"
  compRodSize["dbsnp"]= (`grep -c VC=INDEL #{dbSNPvcf}`).strip.to_i #28833350 
  evalFiles.each{|ff|
    system("java -jar #{gatkPath}/GenomeAnalysisTK.jar -T SelectVariants -R #{refFasta} --variant #{ff} -o #{outputPath}/#{ff.split("\/")[-1]}.indel -selectType INDEL")
  }
  evalFiles.map! {|ff| "#{outputPath}/#{ff.split("\/")[-1]}.indel"}
end

combRodFile=(["-V:"]*evalRods.size).zip(evalRods.zip(evalFiles).map {|rodFile| rodFile.join(" ")}).map {|evalString| evalString.join("")}

#combRodFile2=(["--eval:"]*evalRods.size).zip(evalRods.zip(evalFiles).map {|rodFile| rodFile.join(" ")}).map {|evalString| evalString.join("")}

system("java -jar #{gatkPath}/GenomeAnalysisTK.jar -R #{refFasta} -T VariantsMerger  #{combRodFile.join(" ")} -priority #{evalRods.join(",")} -o #{outputPath}/combined.#{evalRods.join("_")}.#{variantClass}.vcf -setKey set")
system("java -jar #{gatkPath}/GenomeAnalysisTK.jar -T VariantEval -R #{refFasta} -D #{dbSNPvcf} #{selectString.join(" ")} -o #{outputPath}/combeval.#{evalRods.join("_")}.#{variantClass}.report -eval #{outputPath}/combined.#{evalRods.join("_")}.#{variantClass}.vcf --evalModule GenotypeConcordance -l INFO")

#system("java -jar #{gatkPath}/GenomeAnalysisTK.jar  -T VariantEval -R /mnt/scratch0/public/genome/human/hg19.major/hg19.major.fa -D /mnt/scratch0/public/genome/human/hg19/dbsnp/dbsnp_132.vcf #{combRodFile2.join(" ")} -o weval.#{evalRods.join("_")}.#{variantClass}.report -l INFO")

ifid=File.open("#{outputPath}/combeval.#{evalRods.join("_")}.#{variantClass}.report","r")


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



if plotVennFlag
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
    
    weights[comp][combArr.size/2-1] -= sum(knownVector)

    weightByNovelty=[]
    weightAll=[]
    
    weights[comp].each{|w| weightByNovelty.push((w.to_f/scale).round(digitNo))}

    for ii in 1...combArr.size
      weightByNoveltyString+="\""+combArr[ii].join("")+"\" = "+"#{Math.sqrt(weightByNovelty[ii-1])}"
      
      if ii < combArr.size/2
        weightAll.push(weightByNovelty[ii-1]+weightByNovelty[ii-1+combArr.size/2])
        weightAllString+="\""+combArr[ii][1..-1].join("")+"\" = "+"#{weightByNovelty[ii-1]+weightByNovelty[ii-1+combArr.size/2]}"  ##"#{Math.sqrt(weightByNovelty[ii-1]+weightByNovelty[ii-1+combArr.size/2])}"
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
      rInstance.eval_R("plot(Venn(SetNames = c(\"#{evalRods.join("\",\"")}\"),#{weightAllString}),doWeights=TRUE,type=\"circles\",show=list(FaceText=\"weight\",SetLabels=TRUE))") #FaceText=\"weight\"
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
  $stderr.puts("\tall\tknown\tnovel")
  rods.each{|kk,vv|
    if kk[0] == "dbsnp"
      $stderr.puts("#{kk[1]}\t#{vv["all"][ind]}\t#{vv["known"][ind]}\t#{vv["novel"][ind]}")
    end
  }
  ind+=1
  $stderr.puts("")
}

$stderr.puts "#{Time.now} DONE"
exit(0)
