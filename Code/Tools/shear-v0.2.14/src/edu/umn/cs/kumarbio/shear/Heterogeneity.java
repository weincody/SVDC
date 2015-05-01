
package edu.umn.cs.kumarbio.shear;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import java.io.PrintStream;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

/**
 * Walks over the input data set, calculating the number of reads seen for
 * diagnostic purposes.
 * <p>
 * Can also count the number of reads matching a given criterion using read
 * filters (see the --read-filter command line argument). Simplest example of a
 * read-backed analysis.
 * <h2>Input</h2>
 * <p>
 * One or more BAM files.
 * </p>
 * <h2>Output</h2>
 * <p>
 * Number of reads seen.
 * </p>
 * <h2>Examples</h2>
 * 
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CountReads \
 *   -o output.txt \
 *   -I input.bam \
 *   [-L input.intervals]
 * </pre>
 */
@DocumentedGATKFeature(groupName = "Quality Control and Simple Analysis Tools", extraDocs = { CommandLineGATK.class })
@Requires({ DataSource.READS, DataSource.REFERENCE })
public class Heterogeneity extends ReadWalker<Pair<Integer, Integer>, Pair<Integer, Integer>> {

  @Output
  protected PrintStream out;

  @Argument(fullName = "clipLocus", shortName = "clipLocus", doc = "Locus of the breakpoint.", required = true)
  protected int         LOCUS                   = -1;

  @Argument(fullName = "clipDirection", shortName = "clipDirection", doc = "Clip direction from the given locus. Can be + or - for forward or reverse.", required = true)
  protected String      DIRECTION               = "+";

  @Argument(fullName = "unclippedSpanningOnly", shortName = "unclippedSpanningOnly", doc = "Only count spanning reads that are unclipped.", required = false)
  protected boolean     UNCLIPPED_SPANNING_ONLY = false;

  @Override
  public Pair<Integer, Integer> reduceInit() {
    return new Pair<Integer, Integer>(0, 0);
  }

  @Override
  public Pair<Integer, Integer> reduce(Pair<Integer, Integer> value, Pair<Integer, Integer> sum) {
    sum.set(sum.getFirst() + value.getFirst(), sum.getSecond() + value.getSecond());
    return sum;
  }

  @Override
  public Pair<Integer, Integer> map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {

    if (read.getReadUnmappedFlag()) {
      return new Pair<Integer, Integer>(0, 0);
    }

    if (DIRECTION.equals("+")) {
      int pos = read.getAlignmentStart();
      boolean isClipped = false;
      for (CigarElement elt : read.getCigar().getCigarElements()) {
        switch (elt.getOperator()) {
          case H: // hard-clip: ignore
            isClipped = true;
            break;
          case S: // soft clip: return true if we are at the breakpoint
            isClipped = true;
            if (pos == (LOCUS + 1)) {
              // System.out.println(read.getReadName());
              return new Pair<Integer, Integer>(1, 0);
            }
            break;
          case EQ: // sequence match: pass to alignment match
          case X: // sequence mismatch: pass to alignment match
          case M: // alignment match: add to pos
            pos += elt.getLength();
            break;
          case N: // reference skip: pass to deletion
          case D: // deletion: add to pos
            pos += elt.getLength();
            break;
          case P: // padding: don't add to pos
            break;
          case I: // insertion: don't add to pos
            break;
          default:
            throw new IllegalStateException("Case statement didn't deal with cigar op: " + elt.getOperator());
        }
      }
      if (pos == (LOCUS + 1)) {
        return new Pair<Integer, Integer>(0, 0);
      } else if (UNCLIPPED_SPANNING_ONLY && isClipped) {
        return new Pair<Integer, Integer>(0, 0);
      } else {
        return new Pair<Integer, Integer>(0, 1);
      }
    } else if (DIRECTION.equals("-")) {
      if (read.getAlignmentStart() == LOCUS) {
        if (read.getCigar().getCigarElement(0).getOperator() == CigarOperator.S) {
          // System.out.println(read.getReadName());
          return new Pair<Integer, Integer>(1, 0);
        } else {
          return new Pair<Integer, Integer>(0, 0);
        }
      } else {
        boolean isClipped = false;
        for (CigarElement elt : read.getCigar().getCigarElements()) {
          if ((elt.getOperator() == CigarOperator.H) || (elt.getOperator() == CigarOperator.S)) {
            isClipped = true;
          }
        }
        if (UNCLIPPED_SPANNING_ONLY && isClipped) {
          return new Pair<Integer, Integer>(0, 0);
        } else {
          return new Pair<Integer, Integer>(0, 1);
        }
      }
    } else {
      throw new UserException.BadArgumentValue("clipDirection", DIRECTION);
    }
  }

  @Override
  public void onTraversalDone(Pair<Integer, Integer> result) {
    out.print(result.getFirst() + "\n" + result.getSecond() + "\n");
    System.out.println(String.format("Number of soft-clips at the breakpoint: %d", result.getFirst()));
    System.out.println(String.format("Number of reads spanning breakpoint:    %d", result.getSecond()));
  }
}
