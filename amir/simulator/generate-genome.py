import argparse
import random

parser = argparse.ArgumentParser(description='Generate random reference genome.')
parser.add_argument("-c", "--chromosomes",
                  help="Number of chromosomes to generate", metavar="INT", required=True)

parser.add_argument("-l", "--length",
                  help="Chromosome length (each chromosomes lenght will be randomly varied +- length/4)", metavar="INT", required=True)

parser.add_argument("-o", "--output", help="path to output file", metavar="FILE", required=True);

args = parser.parse_args()

out = open(args.output, 'w')

for i in range (0, int(args.chromosomes)):
    chr_name = "binaChr" + str(i)
    out.write(">" + chr_name + "\n")
    random_offset = random.randint(-int(args.length) / 4, int(args.length) / 4)
    length = int(args.length) + random_offset
    bases =  ["A","C", "T", "G"] 
    print (chr_name + ":" + str(length) + "bp")
    for i in range(0, length): 
        #generate random bases
        out.write(bases[random.randint(0,len(bases)-1)])
    out.write("\n")
out.close();

print ("Successfully generated: " + args.output)

