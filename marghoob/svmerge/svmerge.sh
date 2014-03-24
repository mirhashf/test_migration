#!/bin/bash

#cat $* > all.bed

#bedtools sort -i all.bed > sorted.bed

#bedtools merge -i sorted.bed > merged.bed

mkdir -p merged1

merge_opts="-n -d 0 -nms"

for i in raw/*.bed
do 
	echo "First level merging: $i" 1>&2
	merged1="merged1/`basename $i`"
	bedtools sort -i $i | bedtools merge $merge_opts > $merged1;
	bedtools intersect -a $merged1 -b $i -wo > $merged1.details
done
cat raw/*.bed > merged1/.all.bed

mkdir -p merged2

echo "Creating non-overlapping bins"
all=merged2/.all.bed
bins=merged2/.bins.bed
cat merged1/*.bed | bedtools sort > $all
bedtools merge $merge_opts -i $all > $bins

for i in merged1/*.bed
do
	echo "Second level merging: $i" 1>&2
	merged2="merged2"/`basename $i`
	bedtools intersect -f 0.5 -r -a $i -b $bins -wa -u > ${merged2/.bed/.f50r.bed}
	bedtools intersect -f 0.5 -r -a $i -b $bins -v > ${merged2/.bed/.vf50r.bed}
done

cat merged2/*.f50r.bed | bedtools sort | bedtools merge $merge_opts > merged2/.merged.f50r.bed
cat merged2/*.vf50r.bed | bedtools sort | bedtools merge $merge_opts > merged2/.merged.vf50r.bed

bins1=merged2/.merged.f50r.1.bed
bins2=merged2/.merged.f50r.2.bed
awk '{if ($5>1) print $0 > "'$bins2'"; else print $0 > "'$bins1'";}' merged2/.merged.f50r.bed
