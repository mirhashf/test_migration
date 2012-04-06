#!/usr/bin/perl
use strict;

open SNAP_SAM, $ARGV[0] or die $!;
my @prevFields;

while (<SNAP_SAM>) {
  if (m/^@/) {
    print;
    next;
  }

  chomp;
  my @fields = split;
  my $cigarString = $fields[5];  

  my @nameFields = split(/\//, $fields[0]);
  $fields[0] = $nameFields[0];

  if (!($fields[1] & 0x4)) {
    # if the segment is mapped
    my (undef,@cigarOps) = split(/\d+/, $cigarString);
    my @cigarLens = split(/[IDSMX=]/, $cigarString);

    my $transformedCigar = "";
    my $matchLen = 0;
    for (my $cigPos = 0; $cigPos < @cigarOps; $cigPos++) {
      if (($cigarOps[$cigPos] eq "=") or ($cigarOps[$cigPos] eq "X")) {
        $matchLen = $matchLen + $cigarLens[$cigPos];
        next;
      }

      $transformedCigar .= ("$matchLen" . "M") if ($matchLen > 0);
      $transformedCigar .= $cigarLens[$cigPos] . $cigarOps[$cigPos];
      $matchLen = 0;
    }
    $transformedCigar .= ("$matchLen" . "M") if ($matchLen > 0);
    $fields[5] = $transformedCigar;

    # reverse complement, if necessary
#    if ($fields[1] & 0x10) {
#      my $bases = reverse $fields[9];
#      my $quals = reverse $fields[10];

#      $bases =~ tr/ACGTacgt/TGCAtgca/;

#      $fields[9] = $bases;
#      $fields[10] = $quals;
#    }
  } else {
    # if its not mapped, get rid of the reverse complement
    $fields[1] ^= 0x10 if ($fields[1] & 0x10);
  }

  if (@prevFields > 0 and $prevFields[0] eq $fields[0]) {
    $fields[1] ^= 0x8 if (!($prevFields[1] & 0x4) ^ !($fields[1] & 0x8));
    $prevFields[1] ^= 0x8 if (!($fields[1] & 0x4) ^ !($prevFields[1] & 0x8));
    $fields[1] ^= 0x20 if (!($prevFields[1] & 0x10) ^ !($fields[1] & 0x20));
    $prevFields[1] ^= 0x20 if (!($fields[1] & 0x10) ^ !($prevFields[1] & 0x20));

    print join("\t",@prevFields)."\n";
    print join("\t",@fields)."\n";
    undef @prevFields;
  } else {
    print join("\t",@prevFields)."\n" if (@prevFields);
    @prevFields = @fields;
  }
}

print join("\t",@prevFields)."\n" if (@prevFields);
