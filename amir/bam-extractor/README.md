Reads an input BAM file and extracts reads based on specified criteria.

Author: Amirhossein Kiani (amir@binatechnologies.com)

#Building

Run <code> $ ant </code>.

#Running

<pre><code>

$ java -jar target/BamExtractor.jar 
 --in_bam_file <bam_file>  : BAM file to extract reads from.
 --max_length N            : Reads with trimmed-lenght of <= this length will
                             be written to the output file. Default: 150
 --min_length N            : Reads with trimmed-lenght of >= this length will
                             be written to the output file. Defualt: 0
 --out_bam_file <bam_file> : BAM file to write reads to.
 --remove-flags            : Remove RG, PL, PU, LB, SM flags. Default: false
 --trim-quality N          : Quality with which to trim reads. Default: 30

</code></pre>
