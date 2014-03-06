#include <samtools.h>
#include <bam.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SAMPLES 200

int rg_to_int(const char ** rgs, int n_rgs, const char * rg) {
  for (int i = 0; i < n_rgs; i++) {
    if (strcmp(rgs[i], rg) == 0) return i;
  }
  return -1;
}

bamFile init_rg_bam(const char* rg, const char* prefix, const bam_header_t* header) {
  char name[512];
  sprintf(name, "%s.%s.rg.bam", prefix, rg);
  fprintf(stderr, "Opening %s for writing\n", name);
  bamFile out_bam = bam_open(name, "w1");
  bam_header_write(out_bam, header);
  return out_bam;
}

int main(int argc, const char** argv) {
  if (argc < 3) {
    fprintf(stderr, "split_bam input_bam rg1 rg2 ... output_prefix\n");
    exit(1);
  }
  int nsamples = argc - 3;
  int discover = nsamples == 0;

  if (discover) {
    fprintf(stderr, "Will collect rg list from read names\n");
  }

  bamFile inputFile = bam_open(argv[1], "r");
  bam_header_t* header = bam_header_read(inputFile);

  bamFile bam_output_files[MAX_SAMPLES];
  long long counts[MAX_SAMPLES];
  const char * rgs[MAX_SAMPLES];

  const char * output_prefix = argv[argc-1];
  for (int i = 0; i < nsamples; i++) {
    rgs[i] = argv[2+i];
    counts[i] = 0;
    bam_output_files[i] = init_rg_bam(rgs[i], output_prefix, header);
  }

  bam1_t* aln = bam_init1();
  long long total_count = 0;
  while (bam_read1(inputFile, aln) > 0) {
    const char* rg = bam_aux_get(aln, "RG");
    int rg_int = rg_to_int(rgs, nsamples, rg + 1);
    total_count++;
    if (rg_int == -1 && discover) {
      bam_output_files[nsamples] = init_rg_bam(rg + 1, output_prefix, header);
      counts[nsamples] = 0;
      rgs[nsamples] = (char *) malloc(strlen(rg+1) + 1);
      strcpy(rgs[nsamples], rg+1);
      rg_int = nsamples;
      nsamples++;
    } else if (rg_int == -1) {
      continue;
    }
    bam_write1(bam_output_files[rg_int], aln);
    counts[rg_int]++;
  }

  for (int i = 0; i < nsamples; i++) {
    bam_close(bam_output_files[i]);
    printf("%s: %lld reads\n", rgs[i], counts[i]);
  }
  printf("Total: %lld reads\n", total_count);
  return 0;
}
