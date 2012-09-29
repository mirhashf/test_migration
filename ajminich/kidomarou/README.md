h1. Project Kidomarou
by AJ Minich
September-October 2012

h2. Overview

This project seeks to properly encapsulate the ways by 
which one set of genomics files can be processed or 
converted into another set. The goal of the project is 
to simply the way in which processing pipelines are 
understood, configured, and documented. Ultimately, 
this work will become a new implementation of the 
whole-genome sequencing pipeline currently used in the 
Bina Box.

This project is named after Kidomarou, one of the members
of Orochimaru's Sound Four (or Five) in the Naruto anime
series. Kidomarou had spider-like jutsu that allowed him 
to summon lesser spiders, build webs, and see his 
opponents with multiple eyes. His combat strategy involved 
studying his enemies carefully before choosing the optimal 
course of action. The connection here, of course, is that 
this system will navigate through a web of possible 
operations to simplify and optimize the processing of 
genomics data.  

h2. Interface

A user chooses two items:
1. A set of input data to process.
2. A set of desired output data or results.

The system then determines which pipelines are available 
that will process the input data and deliver the desired 
results. There are three cases for these pipelines:
1. There is exactly one pipeline that performs the desired
processing (for example, performing VQSR on variant calls). 
In this case, the user can immediately proceed with the 
only available pipeline.
2. There are several pipelines that perform the desired 
processing, with certain tradeoffs in speed and accuracy.
For example, there are multiple ways to realign and genotype
BAM alignment files. The user is presented with these 
choices of pipelines, and may choose the desired route.
3. There are no pipelines that can perform the desired 
function (for example, turning VCFs back into FASTQ reads).  

h2. Architecture

The system works by using standard graph techniques. Each 
possible set of data files is a node in the graph, 
encapsulated by a FileGroup object or a derived class thereof.
Each possible operation is an edge between these data files, 
with an associated weight depending on how fast the operation 
is. When the user selects an input file group and a desired 
output, the system finds the available paths using typical 
graph pathfinding.

h2. To Do
- install NetworkX
- determine whether there are Python enums
  -> FileDivision (None, Coordinate, Lane)
  -> FileType (FASTA, SAM, BAM, VCF, GFF)
  -> SortOrder (None, Coordinate, Queryname)

h3. Encapsulation

There are still several encapsulation details to work out.
The main issue is that there are two levels of data here:
1. the available nodes and edges in the graph
2. the actual data that can be pushed across the graph

I'm currently having trouble determining how to 
encapsulate these two levels of data in a meaningful, 
easy-to-understand, yet efficient way.
  
h3. Multiple Results

Standard graph techniques are great when there's a 
point-to-point operation, but what about when the user wants 
multiple results (VCFs, GFFs, and AlignStats on some files)?
I haven't figured out how that would be calculated, or how 
those options would be made available to the user.