package io.bina.tools.experimental;

import java.io.File;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Simple utility to compare to reference files
 *
 * @author amir
 */
public class CompareFasta {
  private final static Logger log = LoggerFactory.getLogger(CompareFasta.class);
  @Option(name = "--left", usage = "Left FASTA.", required = true)
  private File left;

  @Option(name = "--right", usage = "Right FASTA.", required = true)
  private File right;

  public static void main(final String[] args) throws Exception {
    new CompareFasta().doMain(args);
  }

  private void doMain(final String[] args) throws Exception {
    final CmdLineParser parser = new CmdLineParser(this);
    parser.setUsageWidth(80);
    try {
      parser.parseArgument(args);

      System.out.println(String.format("chr\tpos\t%s\t%s", left.getName(), right.getName()));
      IndexedFastaSequenceFile leftFasta = new IndexedFastaSequenceFile(left);
      IndexedFastaSequenceFile rightFasta = new IndexedFastaSequenceFile(right);

      ReferenceSequence leftContig;
      boolean theSame = true;
      while ((leftContig = leftFasta.nextSequence()) != null) {
        ReferenceSequence rightContig = rightFasta.nextSequence();
        log.info(String.format("Analyzing contig: %s", leftContig.getName()));

        if (!leftContig.getName().equals(rightContig.getName())) {
          log.error(String.format("Left contig (%s) not eqaual to right contig (%s)", leftContig.getName(), rightContig.getName()));
          break;
        }

        if (leftContig.length() != rightContig.length()) {
          log.error(String.format("Left contig (%s[%i])  not the same size as right contig (%s[%i])", leftContig.getName(), leftContig.length(), rightContig.getName(), rightContig.length()));
          break;
        }
        int pos = 0;
        while (pos < leftContig.getBases().length) {
          if (leftContig.getBases()[pos] != rightContig.getBases()[pos]) {
            theSame = false;
            System.out.println(String.format("%s\t%d\t%s\t%s", leftContig.getName(), pos, (char)leftContig.getBases()[pos], (char)rightContig.getBases()[pos]));
          }
          pos++;
        }
      }
      if (theSame) {
        log.info("The two reference are the same.");
      }

    } catch (final CmdLineException e) {
      parser.printUsage(System.err);
    }
  }
}
