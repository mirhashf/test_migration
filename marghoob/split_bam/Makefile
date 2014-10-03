MMGIT=/net/kodiak/volumes/river/shared/users/marghoob/git
INCLUDES=-I$(MMGIT)/samtools -I$(MMGIT)/htslib
LIBS=-L$(MMGIT)/samtools -L$(MMGIT)/htslib -lbam -lhts

all: split_bam.c
	gcc $(INCLUDES) $(LIBS) split_bam.c -o split_bam -std=c99 -lm -O3 -pthread

clean:
	rm -f split_bam
