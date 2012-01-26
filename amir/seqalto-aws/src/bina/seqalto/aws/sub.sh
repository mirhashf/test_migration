for lane in 1 2 3 4 5 6 7 8
do
	./reads_to_bam.py A804NLABXX.s_$lane seqalto snyder 804 True True m2.4xlarge 50 &
done
