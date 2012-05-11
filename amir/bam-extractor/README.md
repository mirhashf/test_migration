Reads an input BAM file and extracts reads based on specified criteria.

Author: Amirhossein Kiani (amir@binatechnologies.com)
<pre><code>--in_bam_file <bam_file>  : BAM file to extract reads from.
 --max_length N            : Reads with trimmed-lenght of <= this length will
                             be written to the output file. Default: 150
 --min_length N            : Reads with trimmed-lenght of >= this length will
                             be written to the output file. Defualt: 0
 --out_bam_file <bam_file> : BAM file to write reads to.
 --trim-quality N          : Quality with which to trim reads. Default: 30
</code></pre>
