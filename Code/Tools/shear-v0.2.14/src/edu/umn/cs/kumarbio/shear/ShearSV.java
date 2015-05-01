
package edu.umn.cs.kumarbio.shear;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import java.util.TreeSet;
import java.util.concurrent.ExecutionException;

import edu.umn.cs.kumarbio.Genome;
import edu.umn.cs.kumarbio.UnexpectedBehaviorException;

public class ShearSV {

  private String  prefix;
  private String  origBamFile;
  private String  origFaRefFile;
  private String  origBwaIndexPrefix;
  private String  origTwoBitRefFile;
  private String  origPreds;
  private String  region;
  private String  crestDir;
  private String  picardDir;
  private double  minHetThresh;
  private boolean svOnly;
  private boolean hetOnly;
  private boolean debug;
  private String  rgId;
  private String  rgLb;
  private String  rgPl;
  private String  rgPu;
  private String  rgSm;

  /**
   * Creates a new ShearSV object by initializing the local variables that are
   * read from the command line arguments.
   * 
   * @param args
   *          The command line arguments to process, excluding the first one
   *          (SHEAR command, i.e. 'sv' or 'assemble')
   */
  public ShearSV(String[] args) throws ShearUsageException, IOException {

    // Read in arguments
    debug = false;
    String minHetThreshString = null;
    minHetThresh = 0d;
    crestDir = System.getenv("CREST_DIRECTORY");
    picardDir = System.getenv("PICARD_DIRECTORY");
    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-p")) {
        prefix = args[++i];
      } else if (args[i].equals("-b")) {
        origBamFile = args[++i];
      } else if (args[i].equals("-f")) {
        origFaRefFile = args[++i];
      } else if (args[i].equals("-i")) {
        origBwaIndexPrefix = args[++i];
      } else if (args[i].equals("-t")) {
        origTwoBitRefFile = args[++i];
      } else if (args[i].equals("-r")) {
        region = args[++i];
      } else if (args[i].equals("-d")) {
        debug = true;
      } else if (args[i].equals("--preds")) {
        origPreds = args[++i];
      } else if (args[i].equals("--sv-only")) {
        svOnly = true;
      } else if (args[i].equals("--het-only")) {
        hetOnly = true;
      } else if (args[i].equals("--min-het")) {
        minHetThreshString = args[++i];
      } else {
        throw new ShearUsageException("Unknown argument: " + args[i]);
      }
    }

    // Validate prefix
    if (prefix == null) {
      throw new ShearUsageException("Must provide a prefix with the '-p' argument");
    }
    if (!prefix.matches("^[a-zA-Z0-9_\\-\\.]+$")) {
      throw new ShearUsageException("Prefix must consist of alphanumeric characters, underscores, dashes, or periods.");
    }

    // Validate origBamFile
    if (origBamFile == null) {
      throw new ShearUsageException("Must provide a BAM alignment file with the '-b' argument");
    }
    if (!new File(origBamFile).exists()) {
      throw new ShearUsageException("BAM alignment file does not exist: " + origBamFile);
    } else if (!new File(origBamFile + ".bai").exists()) {
      throw new ShearUsageException("BAM alignment index file does not exist: " + origBamFile + ".bai");
    } else {
      origBamFile = new File(origBamFile).getCanonicalPath();
    }

    // Validate origFaRefFile
    if (origFaRefFile == null) {
      throw new ShearUsageException("Must provide a FASTA reference file with the '-f' argument");
    }
    if (!new File(origFaRefFile).exists()) {
      throw new ShearUsageException("FASTA reference file does not exist: " + origFaRefFile);
    } else if (!new File(origFaRefFile + ".fai").exists()) {
      throw new ShearUsageException("FASTA reference index file does not exist: " + origFaRefFile + ".fai");
    } else {
      origFaRefFile = new File(origFaRefFile).getCanonicalPath();
    }

    // Validate origBwaIndexPrefix
    if (origBwaIndexPrefix == null) {
      origBwaIndexPrefix = origFaRefFile;
    }
    if (!new File(origBwaIndexPrefix + ".bwt").exists()) {
      throw new ShearUsageException("BWA index prefix appears to be invalid: " + origBwaIndexPrefix);
    } else {
      origBwaIndexPrefix = new File(origBwaIndexPrefix + ".bwt").getCanonicalPath();
      origBwaIndexPrefix = origBwaIndexPrefix.substring(0, origBwaIndexPrefix.length() - 4);
    }

    // Validate origTwoBitRefFile
    if (origTwoBitRefFile == null) {
      throw new ShearUsageException("Must provide a Two-Bit reference file with the '-t' argument");
    }
    if (!new File(origTwoBitRefFile).exists()) {
      throw new ShearUsageException("Two-Bit reference file does not exist: " + origTwoBitRefFile);
    } else {
      origTwoBitRefFile = new File(origTwoBitRefFile).getCanonicalPath();
    }

    // Validate region
    if (region != null && !region.matches("[^:]+(:[0-9,]+\\-[0-9,]+)?")) {
      throw new ShearUsageException("Invalid region specified: " + region);
    }

    // Validate origPreds
    if (origPreds != null) {
      if (!new File(origPreds).exists()) {
        throw new ShearUsageException("CREST prediction file does not exist: " + origPreds);
      } else {
        origPreds = new File(origPreds).getCanonicalPath();
      }
    }

    // Validate minHetThresh
    if (minHetThreshString != null) {
      try {
        minHetThresh = Double.parseDouble(minHetThreshString);
        if (minHetThresh < 0 || minHetThresh > 1) {
          throw new ShearUsageException("Invalid minimum heterogeneity level specified: " + minHetThreshString);
        }
        minHetThresh = 100d * minHetThresh;
      } catch (NumberFormatException e) {
        throw new ShearUsageException("Invalid minimum heterogeneity level specified: " + minHetThreshString);
      }
    }

    // Validate crestDir
    if (crestDir == null) {
      throw new ShearUsageException("Must specify the CREST directory via the CREST_DIRECTORY environment variable");
    }
    if (!new File(crestDir).exists()) {
      throw new ShearUsageException("CREST directory does not exist: " + crestDir);
    } else {
      crestDir = new File(crestDir).getCanonicalPath();
    }
    if (!new File(crestDir + "/extractSClip.pl").exists()) {
      throw new ShearUsageException("CREST directory does not contain the necessary script extractSClip.pl: "
          + crestDir);
    }
    if (!new File(crestDir + "/CREST.pl").exists()) {
      throw new ShearUsageException("CREST directory does not contain the necessary script CREST.pl: " + crestDir);
    }

    // Validate picardDir
    if (picardDir == null) {
      throw new ShearUsageException(
          "Must specify the Picard Tools directory via the PICARD_DIRECTORY environment variable");
    }
    if (!new File(picardDir).exists()) {
      throw new ShearUsageException("Picard Tools directory does not exist: " + picardDir);
    } else {
      picardDir = new File(picardDir).getCanonicalPath();
    }
    if (!(new File(picardDir + "/picard.jar").exists() || new File(picardDir + "/CleanSam.jar").exists()
        && new File(picardDir + "/FixMateInformation.jar").exists()
        && new File(picardDir + "/AddOrReplaceReadGroups.jar").exists())) {
      throw new ShearUsageException("Picard Tools directory (" + picardDir
          + ") does not contain the necessary JAR file(s), " + "must contain either picard.jar, or all of "
          + "CleanSam.jar, FixMateInformation.jar, and AddOrReplaceReadGroups.jar");
    }

    // Output settings information
    System.out.println("[SHEAR] Starting SHEAR-SV");
    System.out.println("[SHEAR] Prefix:                 " + prefix);
    System.out.println("[SHEAR] BAM alignment file:     " + origBamFile);
    System.out.println("[SHEAR] FASTA reference file:   " + origFaRefFile);
    System.out.println("[SHEAR] BWA index prefix:       " + origBwaIndexPrefix);
    System.out.println("[SHEAR] Two-Bit reference file: " + origTwoBitRefFile);
    System.out.println("[SHEAR] CREST directory:        " + crestDir);
    System.out.println("[SHEAR] Picard Tools directory: " + picardDir);
    if (region != null) {
      System.out.println("[SHEAR] Region:                 " + region);
    }
    if (origPreds != null) {
      System.out.println("[SHEAR] CREST predictions:      " + origPreds);
    }
    if (svOnly) {
      System.out.println("[SHEAR] Prediction mode:        SV Only");
    } else {
      System.out.println("[SHEAR] Prediction mode:        All");
    }
    if (minHetThresh > 0) {
      System.out.printf("[SHEAR] Minimum heterogeneity:  %.2f%%\n", minHetThresh);
    } else {
      System.out.println("[SHEAR] Minimum heterogeneity:  N/A");
    }
    if (debug) {
      System.out.println("[SHEAR] Debugging mode:         ON");
    } else {
      System.out.println("[SHEAR] Debugging mode:         OFF");
    }
    if (hetOnly) {
      System.out.println("[SHEAR] Heterogeneity only:     ON");
    }

  }

  private ArrayList<CrestPrediction> getCrestPreds(String crestPredFile) throws FileNotFoundException {
    ArrayList<CrestPrediction> preds = new ArrayList<CrestPrediction>();
    Scanner sc = new Scanner(new File(crestPredFile));
    while (sc.hasNextLine()) {
      String predLine = sc.nextLine();
      String[] tokens = predLine.split("\t");
      // What about cross-chromosome events?
      String bp1region = tokens[0];
      int bp1 = Integer.parseInt(tokens[1]);
      String bp1dir = tokens[2];
      String bp2region = tokens[4];
      int bp2 = Integer.parseInt(tokens[5]);
      String bp2dir = tokens[6];
      String type = tokens[8];
      preds.add(new CrestPrediction(type, bp1region, bp1, bp1dir, 0, 0, bp2region, bp2, bp2dir, 0, 0, 0));
    }
    sc.close();
    return preds;
  }

  private ArrayList<CrestPrediction> filterCrestPreds(ArrayList<CrestPrediction> origPreds) {
    // SHOULD MAKE SURE SORTED
    // NEED TO CHECK FOR SVs WITH SAME LOC ON BOTH BREAKPOINTS
    ArrayList<CrestPrediction> newPreds = new ArrayList<CrestPrediction>();
    tryNextPred: for (int i = 0; i < origPreds.size(); i++) {
      CrestPrediction origPred = origPreds.get(i);
      for (int k = 0; k < newPreds.size(); k++) {
        CrestPrediction newPred = newPreds.get(k);
        if (origPred.bp1 == newPred.bp1 && origPred.bp1region.equals(newPred.bp1region)) {
          if (Math.abs(origPred.bp2 - newPred.bp2) <= 10) {
            if (newPred.bp2sc < origPred.bp2sc) {
              newPreds.set(k, origPred);
            }
            continue tryNextPred;
          }
        } else if (origPred.bp2 == newPred.bp2 && origPred.bp2region.equals(newPred.bp2region)) {
          if (Math.abs(origPred.bp1 - newPred.bp1) <= 10) {
            if (newPred.bp1sc < origPred.bp1sc) {
              newPreds.set(k, origPred);
            }
            continue tryNextPred;
          }
        }
      }
      newPreds.add(origPreds.get(i));
    }
    return newPreds;
  }

  private void outputResults(TreeSet<Variant> variants) throws IOException {

    Genome gen = null;
    PrintWriter allSdiOut = null;
    PrintWriter allReportOut = null;
    PrintWriter snpIndelSdiOut = null;
    PrintWriter snpIndelReportOut = null;
    PrintWriter svSdiOut = null;
    PrintWriter svReportOut = null;
    boolean foundSvs = false;
    boolean foundSnpIndels = false;

    try {

      String reportHeaderLine =
          String.format(
              "%-8s%-30s%-21s%-21s%13s",
              "Type",
              "Location / Size",
              "Reference Bases",
              "Variant Bases",
              "Heterogeneity");
      gen = new Genome(origFaRefFile);
      if (!svOnly) {
        allSdiOut = new PrintWriter(prefix + ".all.sdi");
        allReportOut = new PrintWriter(prefix + ".all.report");
        snpIndelSdiOut = new PrintWriter(prefix + ".snpindel.sdi");
        snpIndelReportOut = new PrintWriter(prefix + ".snpindel.report");
      }
      svSdiOut = new PrintWriter(prefix + ".sv.sdi");
      svReportOut = new PrintWriter(prefix + ".sv.report");

      for (Variant v : variants) {
        // Fix report format to account for chromosome name lengths
        String sdiLine = null;
        String reportLine = null;

        // Output SNP
        if (v instanceof SNP) {
          SNP snp = (SNP) v;
          sdiLine =
              getSdiLine(
                  snp.getChr(),
                  snp.getLoc(),
                  Character.toString(snp.getRefBase()),
                  Character.toString(snp.getVarBase()));
          reportLine =
              String.format(
                  "%-8s%-30s%-21s%-21s%12.2f%%",
                  "SNP",
                  snp.getChr() + ":" + snp.getLoc(),
                  snp.getRefBase(),
                  snp.getVarBase(),
                  snp.getHeterogeneity());

          // Output Deletion
        } else if (v instanceof Deletion) {
          Deletion del = (Deletion) v;
          String deletedBases = gen.getSequence(del.getChr(), del.getLoc(), del.getEndLoc());
          String truncatedDeletedBases =
              deletedBases.length() > 19 ? deletedBases.substring(0, 8) + "..."
                  + deletedBases.substring(deletedBases.length() - 8) : deletedBases;
          sdiLine = getSdiLine(del.getChr(), del.getLoc(), deletedBases, "-");
          reportLine =
              String.format(
                  "%-8s%-30s%-21s%-21s%12.2f%%",
                  "DEL",
                  del.getChr() + ":" + del.getLoc() + "-" + del.getEndLoc(),
                  truncatedDeletedBases,
                  "-",
                  del.getHeterogeneity());

          // Output Insertion
        } else if (v instanceof Insertion) {
          Insertion ins = (Insertion) v;
          String truncatedInsertedBases =
              ins.getInsertedBases().length() > 19 ? ins.getInsertedBases().substring(0, 8) + "..."
                  + ins.getInsertedBases().substring(ins.getInsertedBases().length() - 8) : ins.getInsertedBases();
          sdiLine = getSdiLine(ins.getChr(), ins.getLoc(), "-", ins.getInsertedBases());
          reportLine =
              String.format(
                  "%-8s%-30s%-21s%-21s%12.2f%%",
                  "INS",
                  ins.getChr() + ":" + ins.getLoc() + " (" + ins.getSize() + " bp)",
                  "-",
                  truncatedInsertedBases,
                  ins.getHeterogeneity());

          // Output TandemDuplication
        } else if (v instanceof TandemDuplication) {
          TandemDuplication tdup = (TandemDuplication) v;
          sdiLine =
              getSdiLine(
                  tdup.getChr(),
                  tdup.getLoc(),
                  "-",
                  gen.getSequence(tdup.getChr(), tdup.getLoc(), tdup.getEndLoc()));
          reportLine =
              String.format(
                  "%-8s%-30s%-21s%-21s%12.2f%%",
                  "TDUP",
                  tdup.getChr() + ":" + tdup.getLoc() + "-" + tdup.getEndLoc(),
                  "",
                  "",
                  tdup.getHeterogeneity());

          // Output Inversion
        } else if (v instanceof Inversion) {
          Inversion inv = (Inversion) v;
          String referenceBases = gen.getSequence(inv.getChr(), inv.getLoc(), inv.getEndLoc());
          String invertedBases = Genome.rc(gen.getSequence(inv.getChr(), inv.getLoc(), inv.getEndLoc()));
          sdiLine = getSdiLine(inv.getChr(), inv.getLoc(), referenceBases, invertedBases);
          String truncatedReferenceBases =
              referenceBases.length() > 19 ? referenceBases.substring(0, 8) + "..."
                  + referenceBases.substring(referenceBases.length() - 8) : referenceBases;
          String truncatedInvertedBases =
              invertedBases.length() > 19 ? invertedBases.substring(0, 8) + "..."
                  + invertedBases.substring(invertedBases.length() - 8) : invertedBases;
          reportLine =
              String.format(
                  "%-8s%-30s%-21s%-21s%12.2f%%",
                  "INV",
                  inv.getChr() + ":" + inv.getLoc() + "-" + inv.getEndLoc(),
                  truncatedReferenceBases,
                  truncatedInvertedBases,
                  inv.getHeterogeneity());

          // Output Translocation
        } else if (v instanceof Translocation) {
          Translocation xloc = (Translocation) v;
          reportLine =
              String.format(
                  "%-8s%-30s%-21s%-21s%12.2f%%",
                  "XLOC",
                  xloc.getChr() + ":" + xloc.getLoc() + "-" + xloc.getEndChr() + ":" + xloc.getEndLoc(),
                  "",
                  "",
                  xloc.getHeterogeneity());

          // Unknown variant type
        } else {
          throw new IllegalArgumentException("Encountered an invalid SV type while outputting results");
        }

        // Output to standard output report
        if (!foundSnpIndels && !foundSvs) {
          System.out.println("[SHEAR]     " + reportHeaderLine);
        }
        System.out.println("[SHEAR]     " + reportLine);

        // Output to combined report (only in 'SV only' mode)
        if (!svOnly) {
          if (!foundSnpIndels && !foundSvs) {
            allReportOut.println(reportHeaderLine);
          }
          allReportOut.println(reportLine);
          if (sdiLine != null) {
            allSdiOut.print(sdiLine);
          }
        }

        // Output to SNPINDEL report
        if (!svOnly && (v instanceof SNP || v.isIndel())) {
          if (!foundSnpIndels) {
            snpIndelReportOut.println(reportHeaderLine);
          }
          snpIndelReportOut.println(reportLine);
          foundSnpIndels = true;
          if (sdiLine != null) {
            snpIndelSdiOut.print(sdiLine);
          }

          // Output to SV report
        } else {
          if (!foundSvs) {
            svReportOut.println(reportHeaderLine);
          }
          svReportOut.println(reportLine);
          foundSvs = true;
          if (sdiLine != null) {
            svSdiOut.print(sdiLine);
          }
        }

      }

      // If no results were output
      if (!foundSnpIndels && !foundSvs) {
        System.out.println("[SHEAR] No results!");
      }

    } finally {
      gen.close();
      if (!svOnly) {
        allSdiOut.close();
        allReportOut.close();
        snpIndelSdiOut.close();
        snpIndelReportOut.close();
      }
      svSdiOut.close();
      svReportOut.close();
    }

  }

  private ArrayList<CrestPrediction> processCrestPreds(
      String bamFileName, String refFileName, ArrayList<CrestPrediction> crestPreds) throws Exception {
    for (int i = 0; i < crestPreds.size(); i++) {
      CrestPrediction pred = crestPreds.get(i);
      System.out.printf(
          "[SHEAR] Processing: %s %s:%d(%s) - %s:%d(%s)\n",
          pred.type,
          pred.bp1region,
          pred.bp1,
          pred.bp1dir,
          pred.bp2region,
          pred.bp2,
          pred.bp2dir);
      if (pred.type.equals("DEL")) {
        crestPreds.set(i, processDel(bamFileName, refFileName, pred));
      } else if (pred.type.equals("INS")) {
        crestPreds.set(i, processIns(bamFileName, refFileName, pred));
      } else if (pred.type.equals("ITX")) {
        crestPreds.set(i, processItx(bamFileName, refFileName, pred));
      } else if (pred.type.equals("CTX")) {
        crestPreds.set(i, processCtx(bamFileName, refFileName, pred));
      } else if (pred.type.equals("INV")) {
        crestPreds.set(i, processInv(bamFileName, refFileName, pred));
      } else {
        throw new IllegalArgumentException("The following SV type is not valid: " + pred.type);
      }
    }
    return crestPreds;
  }

  private CrestPrediction processDel(String bamFileName, String refFileName, CrestPrediction pred) throws Exception {

    Random rand = new Random(System.currentTimeMillis());
    String tmpFile1 = "tmp-" + pred.bp1region + pred.bp1 + "-" + rand.nextInt(Integer.MAX_VALUE);
    String tmpFile2 = "tmp-" + pred.bp2region + pred.bp2 + "-" + rand.nextInt(Integer.MAX_VALUE);

    String[] gatkArgs1 =
        { "-T", "Heterogeneity", "-I", bamFileName, "-R", refFileName, "-clipLocus", String.valueOf(pred.bp1),
            "-clipDirection", "+", "-L", pred.bp1region + ":" + pred.bp1, "-o", tmpFile1, "-l", "ERROR" };
    if (debug) {
      gatkArgs1[gatkArgs1.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs1, debug);

    String[] gatkArgs2 =
        { "-T", "Heterogeneity", "-I", bamFileName, "-R", refFileName, "-clipLocus", String.valueOf(pred.bp2),
            "-clipDirection", "-", "-L", pred.bp2region + ":" + pred.bp2, "-o", tmpFile2, "-l", "ERROR" };
    if (debug) {
      gatkArgs2[gatkArgs2.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs2, debug);

    File bp1file = new File(tmpFile1);
    Scanner bp1scanner = new Scanner(bp1file);
    int bp1sc = bp1scanner.nextInt();
    int bp1span = bp1scanner.nextInt();
    bp1scanner.close();
    bp1file.delete();

    File bp2file = new File(tmpFile2);
    Scanner bp2scanner = new Scanner(bp2file);
    int bp2sc = bp2scanner.nextInt();
    int bp2span = bp2scanner.nextInt();
    bp2scanner.close();
    bp2file.delete();

    double per =
        100 * ((double) bp1sc + (double) bp2sc)
            / ((double) bp1sc + (double) bp2sc + ((double) bp1span + (double) bp2span) / 2);
    if (per > 100) {
      per = 100;
    }
    pred.bp1sc = bp1sc;
    pred.bp1span = bp1span;
    pred.bp2sc = bp2sc;
    pred.bp2span = bp2span;
    pred.heterogeneityPercent = per;

    System.out.println(String.format(
        "[SHEAR]     Estimated heterogeneity:                %.2f%%",
        pred.heterogeneityPercent));

    return pred;

  }

  private CrestPrediction processIns(String bamFileName, String refFileName, CrestPrediction pred) throws Exception {

    Random rand = new Random(System.currentTimeMillis());
    String tmpFile1 = "tmp-" + pred.bp1region + pred.bp1 + "-" + rand.nextInt(Integer.MAX_VALUE);
    String tmpFile2 = "tmp-" + pred.bp2region + pred.bp2 + "-" + rand.nextInt(Integer.MAX_VALUE);

    String[] gatkArgs1 =
        { "-T", "Heterogeneity", "-I", bamFileName, "-R", refFileName, "-clipLocus", String.valueOf(pred.bp2),
            "-clipDirection", "-", "-L", pred.bp2region + ":" + pred.bp2, "-o", tmpFile1, "-l", "ERROR" };
    if (debug) {
      gatkArgs1[gatkArgs1.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs1, debug);

    String[] gatkArgs2 =
        { "-T", "Heterogeneity", "-I", bamFileName, "-R", refFileName, "-clipLocus", String.valueOf(pred.bp1),
            "-clipDirection", "+", "-L", pred.bp1region + ":" + pred.bp1, "-o", tmpFile2, "-l", "ERROR" };
    if (debug) {
      gatkArgs2[gatkArgs2.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs2, debug);

    File bp1file = new File(tmpFile1);
    Scanner bp1scanner = new Scanner(bp1file);
    int bp1sc = bp1scanner.nextInt();
    int bp1span = bp1scanner.nextInt();
    bp1scanner.close();
    bp1file.delete();

    File bp2file = new File(tmpFile2);
    Scanner bp2scanner = new Scanner(bp2file);
    int bp2sc = bp2scanner.nextInt();
    int bp2span = bp2scanner.nextInt();
    bp2scanner.close();
    bp2file.delete();

    double per = 100 * ((double) bp1sc + (double) bp2sc) / (((double) bp1span + (double) bp2span) / 2);
    if (per > 100) {
      per = 100;
    }
    pred.bp1sc = bp1sc;
    pred.bp1span = bp1span;
    pred.bp2sc = bp2sc;
    pred.bp2span = bp2span;
    pred.heterogeneityPercent = per;

    System.out.println(String.format(
        "[SHEAR]     Estimated heterogeneity:                %.2f%%",
        pred.heterogeneityPercent));

    return pred;

  }

  private CrestPrediction processItx(String bamFileName, String refFileName, CrestPrediction pred) throws Exception {

    // FIX clip directions

    Random rand = new Random(System.currentTimeMillis());
    String tmpFile1 = "tmp-" + pred.bp1region + pred.bp1 + "-" + rand.nextInt(Integer.MAX_VALUE);
    String tmpFile2 = "tmp-" + pred.bp2region + pred.bp2 + "-" + rand.nextInt(Integer.MAX_VALUE);

    String[] gatkArgs1 =
        { "-T", "Heterogeneity", "-I", bamFileName, "-R", refFileName, "-clipLocus", String.valueOf(pred.bp1),
            "-clipDirection", "+", "-L", pred.bp1region + ":" + pred.bp1, "-o", tmpFile1, "-l", "ERROR" };
    if (debug) {
      gatkArgs1[gatkArgs1.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs1, debug);

    String[] gatkArgs2 =
        { "-T", "Heterogeneity", "-I", bamFileName, "-R", refFileName, "-clipLocus", String.valueOf(pred.bp2),
            "-clipDirection", "-", "-L", pred.bp2region + ":" + pred.bp2, "-o", tmpFile2, "-l", "ERROR" };
    if (debug) {
      gatkArgs2[gatkArgs2.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs2, debug);

    File bp1file = new File(tmpFile1);
    Scanner bp1scanner = new Scanner(bp1file);
    int bp1sc = bp1scanner.nextInt();
    int bp1span = bp1scanner.nextInt();
    bp1scanner.close();
    bp1file.delete();

    File bp2file = new File(tmpFile2);
    Scanner bp2scanner = new Scanner(bp2file);
    int bp2sc = bp2scanner.nextInt();
    int bp2span = bp2scanner.nextInt();
    bp2scanner.close();
    bp2file.delete();

    double per = 100 * ((double) bp1sc + (double) bp2sc) / ((double) bp1sc + (double) bp2sc + bp1span + bp2span);
    if (per > 100) {
      per = 100;
    }
    pred.bp1sc = bp1sc;
    pred.bp1span = bp1span;
    pred.bp2sc = bp2sc;
    pred.bp2span = bp2span;
    pred.heterogeneityPercent = per;

    System.out.println(String.format(
        "[SHEAR]     Estimated heterogeneity:                %.2f%%",
        pred.heterogeneityPercent));

    return pred;

  }

  private CrestPrediction processCtx(String bamFileName, String refFileName, CrestPrediction pred) throws Exception {

    // FIX clip directions

    Random rand = new Random(System.currentTimeMillis());
    String tmpFile1 = "tmp-" + pred.bp1region + pred.bp1 + "-" + rand.nextInt(Integer.MAX_VALUE);
    String tmpFile2 = "tmp-" + pred.bp2region + pred.bp2 + "-" + rand.nextInt(Integer.MAX_VALUE);

    String[] gatkArgs1 =
        { "-T", "Heterogeneity", "-I", bamFileName, "-R", refFileName, "-clipLocus", String.valueOf(pred.bp1),
            "-clipDirection", "+", "-L", pred.bp1region + ":" + pred.bp1, "-o", tmpFile1, "-l", "ERROR" };
    if (debug) {
      gatkArgs1[gatkArgs1.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs1, debug);

    String[] gatkArgs2 =
        { "-T", "Heterogeneity", "-I", bamFileName, "-R", refFileName, "-clipLocus", String.valueOf(pred.bp2),
            "-clipDirection", "-", "-L", pred.bp2region + ":" + pred.bp2, "-o", tmpFile2, "-l", "ERROR" };
    if (debug) {
      gatkArgs2[gatkArgs2.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs2, debug);

    File bp1file = new File(tmpFile1);
    Scanner bp1scanner = new Scanner(bp1file);
    int bp1sc = bp1scanner.nextInt();
    int bp1span = bp1scanner.nextInt();
    bp1scanner.close();
    bp1file.delete();

    File bp2file = new File(tmpFile2);
    Scanner bp2scanner = new Scanner(bp2file);
    int bp2sc = bp2scanner.nextInt();
    int bp2span = bp2scanner.nextInt();
    bp2scanner.close();
    bp2file.delete();

    double per = 100 * ((double) bp1sc + (double) bp2sc) / ((double) bp1sc + (double) bp2sc + bp1span + bp2span);
    if (per > 100) {
      per = 100;
    }
    pred.bp1sc = bp1sc;
    pred.bp1span = bp1span;
    pred.bp2sc = bp2sc;
    pred.bp2span = bp2span;
    pred.heterogeneityPercent = per;

    System.out.println(String.format(
        "[SHEAR]     Estimated heterogeneity:                %.2f%%",
        pred.heterogeneityPercent));

    return pred;

  }

  private CrestPrediction processInv(String bamFileName, String refFileName, CrestPrediction pred) throws Exception {

    Random rand = new Random(System.currentTimeMillis());
    String tmpFile1 = "tmp-" + pred.bp1region + pred.bp1 + "-" + rand.nextInt(Integer.MAX_VALUE);
    String tmpFile2 = "tmp-" + pred.bp2region + pred.bp2 + "-" + rand.nextInt(Integer.MAX_VALUE);

    String[] gatkArgs1 =
        { "-T", "Heterogeneity", "-I", bamFileName, "-R", refFileName, "-clipLocus", String.valueOf(pred.bp1),
            "-clipDirection", "+", "-unclippedSpanningOnly", "-L", pred.bp1region + ":" + pred.bp1, "-o", tmpFile1,
            "-l", "ERROR" };
    if (debug) {
      gatkArgs1[gatkArgs1.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs1, debug);

    String[] gatkArgs2 =
        { "-T", "Heterogeneity", "-I", bamFileName, "-R", refFileName, "-clipLocus", String.valueOf(pred.bp2),
            "-clipDirection", "-", "-unclippedSpanningOnly", "-L", pred.bp2region + ":" + pred.bp2, "-o", tmpFile2,
            "-l", "ERROR" };
    if (debug) {
      gatkArgs2[gatkArgs2.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs2, debug);

    File bp1file = new File(tmpFile1);
    Scanner bp1scanner = new Scanner(bp1file);
    int bp1sc = bp1scanner.nextInt();
    int bp1span = bp1scanner.nextInt();
    bp1scanner.close();
    bp1file.delete();

    File bp2file = new File(tmpFile2);
    Scanner bp2scanner = new Scanner(bp2file);
    int bp2sc = bp2scanner.nextInt();
    int bp2span = bp2scanner.nextInt();
    bp2scanner.close();
    bp2file.delete();

    double per =
        100 * (2d * ((double) bp1sc + (double) bp2sc)) / (2d * ((double) bp1sc + (double) bp2sc) + bp1span + bp2span);
    if (per > 100) {
      per = 100;
    }
    pred.bp1sc = bp1sc;
    pred.bp1span = bp1span;
    pred.bp2sc = bp2sc;
    pred.bp2span = bp2span;
    pred.heterogeneityPercent = per;

    System.out.println(String.format(
        "[SHEAR]     Estimated heterogeneity:                %.2f%%",
        pred.heterogeneityPercent));

    return pred;

  }

  /*
   * Returns the properly formatted SDI line for a given region name, start
   * location, reference string, and consensus string.
   */
  public static String getSdiLine(String regionName, int startLoc, String reference, String consensus) {
    StringBuffer line = new StringBuffer();
    line.append(regionName);
    line.append("\t");
    line.append(startLoc);
    line.append("\t");
    if (reference.equals("-")) {
      if (consensus.equals("-")) {
        line.append(0);
      } else {
        line.append(consensus.length());
      }
    } else if (consensus.equals("-")) {
      if (reference.equals("-")) {
        line.append(0);
      } else {
        line.append("-");
        line.append(reference.length());
      }
    } else {
      line.append(consensus.length() - reference.length());
    }
    line.append("\t");
    line.append(reference);
    line.append("\t");
    line.append(consensus);
    line.append("\t*\t*\n");
    return line.toString();
  }

  /*
   * Prepare the original BAM file by running CleanSam, FixMateInformation,
   * sorting, and ensuring read groups are listed. Also saves the read groups
   * for future use.
   */
  private void prepareOrigBamFile(String bamFile)
      throws InterruptedException, IOException, ExecutionException, Exception {

    // Temp output string
    String tmpOutput = prefix + ".tmp" + System.currentTimeMillis();

    // Get the headers from the original BAM file
    FileOutputStream output = new FileOutputStream(tmpOutput + ".origHeaders");
    if (debug) {
      RunSystemCommand.run(
          new ProcessBuilder("samtools", "view", "-H", bamFile),
          output,
          null,
          System.out,
          "samtools",
          debug);
    } else {
      RunSystemCommand.run(new ProcessBuilder("samtools", "view", "-H", bamFile), output, null, null, null, debug);
    }
    output.close();

    // Read the headers to find read group tags
    Scanner sc = new Scanner(new File(tmpOutput + ".origHeaders"));
    rgId = prefix;
    rgLb = prefix;
    rgPl = prefix;
    rgPu = prefix;
    rgSm = prefix;
    while (sc.hasNextLine()) {
      String line = sc.nextLine();
      if (line.startsWith("@RG")) {
        String[] tokens = line.split("\t");
        for (int i = 1; i < tokens.length; i++) {
          if (tokens[i].startsWith("ID")) {
            rgId = tokens[i].substring(3);
          } else if (tokens[i].startsWith("LB")) {
            rgLb = tokens[i].substring(3);
          } else if (tokens[i].startsWith("PL")) {
            rgPl = tokens[i].substring(3);
          } else if (tokens[i].startsWith("PU")) {
            rgPu = tokens[i].substring(3);
          } else if (tokens[i].startsWith("SM")) {
            rgSm = tokens[i].substring(3);
          }
        }
        if (rgLb == null) {
          rgLb = rgId;
        }
        if (rgPl == null) {
          rgPl = rgId;
        }
        if (rgPu == null) {
          rgPu = rgId;
        }
        if (rgSm == null) {
          rgSm = rgId;
        }
        break;
      }
    }
    sc.close();

    // Remove temp headers file
    new File(tmpOutput + ".origHeaders").delete();

    // Clean BAM file
    String[] picardArgs1 = { "I=" + bamFile, "O=" + tmpOutput + ".cleaned.bam", "QUIET=TRUE", "VERBOSITY=ERROR" };
    if (debug) {
      picardArgs1[picardArgs1.length - 2] = "QUIET=FALSE";
      picardArgs1[picardArgs1.length - 1] = "VERBOSITY=INFO";
    }
    RunPicardCommand.run(picardDir, "CleanSam", picardArgs1, debug);

    // Fix mate information for BAM file
    String[] picardArgs2 =
        { "I=" + tmpOutput + ".cleaned.bam", "O=" + tmpOutput + ".cleaned.fixed.bam", "SO=coordinate", "QUIET=TRUE",
            "VERBOSITY=ERROR" };
    if (debug) {
      picardArgs2[picardArgs2.length - 2] = "QUIET=FALSE";
      picardArgs2[picardArgs2.length - 1] = "VERBOSITY=INFO";
    }
    RunPicardCommand.run(picardDir, "FixMateInformation", picardArgs2, debug);

    // Add the read group info back, in case it was missing before
    addReadGroupsAndSort(tmpOutput + ".cleaned.fixed.bam", prefix + "-orig-with-rg.bam");
    new File(tmpOutput + ".cleaned.bam").delete();
    new File(tmpOutput + ".cleaned.fixed.bam").delete();

    // Index the prepared BAM file
    if (debug) {
      RunSystemCommand.run(
          new ProcessBuilder("samtools", "index", prefix + "-orig-with-rg.bam"),
          System.out,
          "samtools",
          System.out,
          "samtools",
          debug);
    } else {
      RunSystemCommand.run(
          new ProcessBuilder("samtools", "index", prefix + "-orig-with-rg.bam"),
          null,
          null,
          null,
          null,
          debug);
    }

  }

  private TreeSet<Variant> runHaplotypeCaller(String bamFile, String faRefFile, String vcfOutputFile) throws Exception {
    String[] gatkArgs =
        { "-T", "HaplotypeCaller", "-I", bamFile, "-R", faRefFile, "-o", vcfOutputFile, "-rf", "BadCigar", "-l",
            "ERROR" };
    if (debug) {
      gatkArgs[gatkArgs.length - 1] = "INFO";
    }
    RunGatkCommand.run(gatkArgs, debug);
    Scanner sc = new Scanner(new File(vcfOutputFile));
    TreeSet<Variant> variants = new TreeSet<Variant>();
    while (sc.hasNextLine()) {
      String line = sc.nextLine();
      if (line.charAt(0) != '#') {
        String[] tokens = line.split("\t");
        String chr = tokens[0];
        int loc = Integer.valueOf(tokens[1]);
        String refBases = tokens[3];
        String varBases = tokens[4];
        int variantSize = varBases.length() - refBases.length();
        Variant v = null;
        if (refBases.matches("^[ACGT]+$") && varBases.matches("^[ACGT]+$")) {
          if (variantSize == 0 & refBases.length() == 1) {
            v = new SNP(chr, loc, refBases.charAt(0), varBases.charAt(0));
          } else if (variantSize > 0 & refBases.length() == 1 & refBases.charAt(0) == varBases.charAt(0)) {
            // Check for weird case VCF-locs (ends of starts of chr)
            v = new Insertion(chr, loc + 1, varBases.substring(1));
            v.setIndel(true);
          } else if (variantSize < 0 & varBases.length() == 1 & refBases.charAt(0) == varBases.charAt(0)) {
            // Check for weird case VCF-locs (ends of starts of chr)
            v = new Deletion(chr, loc + 1, -variantSize);
            v.setIndel(true);
          } else {
            if (debug) {
              System.out.println("[SHEAR] Skipping unexpected event in HaplotypeCaller output (" + chr + ":" + loc
                  + " " + refBases + " -> " + varBases + ")");
            }
            continue;
          }
        } else {
          if (debug) {
            System.out.println("[SHEAR] Skipping unexpected event in HaplotypeCaller output (" + chr + ":" + loc + " "
                + refBases + " -> " + varBases + ")");
          }
          continue;
        }
        String[] formatTags = tokens[8].split(":");
        String[] formatInfo = tokens[9].split(":");
        double heterogeneity = Double.NaN;
        for (int i = 0; i < formatTags.length; i++) {
          if (formatTags[i].equals("AD")) {
            String[] alleleDepths = formatInfo[i].split(",");
            double varDepth = Double.valueOf(alleleDepths[1]);
            double refDepth = Double.valueOf(alleleDepths[0]);
            heterogeneity = 100.0 * (varDepth / (varDepth + refDepth));
          }
        }
        if (heterogeneity != Double.NaN) {
          if (heterogeneity > 0) {
            v.setHeterogeneity(heterogeneity);
            variants.add(v);
          } else if (debug) {
            System.out.println("[SHEAR] Filtering out SNP/INDEL (" + chr + ":" + loc + " " + refBases + " -> "
                + varBases + ") with 0% heterogeneity, considered false positive");
          }
        } else {
          sc.close();
          throw new UnexpectedBehaviorException("AD not found or not containing valid data");
        }
      }
    }
    sc.close();
    return variants;
  }

  private void runCrest(
      String bamFileName, String crestPrefix, String faRefFile, GFServer gfServer, boolean onExtractedAlignment)
      throws InterruptedException, IOException, ExecutionException {
    ProcessBuilder extractSClipProcess =
        new ProcessBuilder(
            "perl",
            "-I" + crestDir,
            crestDir + "/extractSClip.pl",
            "-p",
            crestPrefix,
            "-i",
            bamFileName,
            "--ref_genome",
            faRefFile);
    ProcessBuilder crestProcess;
    if (onExtractedAlignment) {
      crestProcess =
          new ProcessBuilder(
              "perl",
              "-I" + crestDir,
              crestDir + "/CREST.pl",
              "--min_sclip_reads",
              "2",
              "--max_rep_cover",
              "1000000",
              "-p",
              crestPrefix,
              "-f",
              crestPrefix + ".cover",
              "-d",
              bamFileName,
              "--ref_genome",
              faRefFile,
              "-t",
              gfServer.getTwoBitFile(),
              "--blatserver",
              "localhost",
              "--blatport",
              "" + gfServer.getPortNum());
    } else {
      crestProcess =
          new ProcessBuilder(
              "perl",
              "-I" + crestDir,
              crestDir + "/CREST.pl",
              "--max_rep_cover",
              "1000000",
              "-p",
              crestPrefix,
              "-f",
              crestPrefix + ".cover",
              "-d",
              bamFileName,
              "--ref_genome",
              faRefFile,
              "-t",
              gfServer.getTwoBitFile(),
              "--blatserver",
              "localhost",
              "--blatport",
              "" + gfServer.getPortNum());
    }
    if (debug) {
      RunSystemCommand.run(extractSClipProcess, System.out, "CREST", System.out, "CREST", debug);
      RunSystemCommand.run(crestProcess, System.out, "CREST", System.out, "CREST", debug);
    } else {
      RunSystemCommand.run(extractSClipProcess, null, null, null, null, debug);
      RunSystemCommand.run(crestProcess, null, null, null, null, debug);
      new File(crestPrefix + ".sclip.txt").delete();
      new File(crestPrefix + ".cover").delete();
    }
  }

  private ArrayList<String> getBreakpoints(String crestPredFile) throws FileNotFoundException {
    Scanner sc = new Scanner(new File(crestPredFile));
    ArrayList<String> breakpoints = new ArrayList<String>();
    while (sc.hasNextLine()) {
      String line = sc.nextLine();
      String[] tokens = line.split("\t");
      breakpoints.add(tokens[0] + ":" + tokens[1]);
      breakpoints.add(tokens[4] + ":" + tokens[5]);
    }
    sc.close();
    return breakpoints;
  }

  public void addReadGroupsAndSort(String inputFile, String outputFile) throws IOException, Exception {
    String[] picardArgs =
        { "I=" + inputFile, "O=" + outputFile, "SO=coordinate", "RGID=" + rgId, "RGLB=" + rgLb, "RGPL=" + rgPl,
            "RGPU=" + rgPu, "RGSM=" + rgSm, "QUIET=TRUE", "VERBOSITY=ERROR" };
    if (debug) {
      picardArgs[picardArgs.length - 2] = "QUIET=FALSE";
      picardArgs[picardArgs.length - 1] = "VERBOSITY=INFO";
    }
    RunPicardCommand.run(picardDir, "AddOrReplaceReadGroups", picardArgs, debug);
  }

  /*
   * private void extractRelevantReads(ArrayList<String> breakpoints, String
   * bamFile, String faRefFile, String fqOutputFile) throws
   * InterruptedException, IOException, Exception { String tmpOutput = prefix +
   * ".tmp" + System.currentTimeMillis(); if (breakpoints.size() > 0) { String[]
   * gatkArgs = new String[12 + (2 * breakpoints.size())]; int ind = 0;
   * gatkArgs[ind++] = "-T"; gatkArgs[ind++] = "PrintReads"; gatkArgs[ind++] =
   * "-I"; gatkArgs[ind++] = bamFile; gatkArgs[ind++] = "-R"; gatkArgs[ind++] =
   * faRefFile; gatkArgs[ind++] = "--interval_padding"; gatkArgs[ind++] = "25";
   * for (int i = 0; i < breakpoints.size(); i++) { gatkArgs[ind++] =
   * "--intervals"; gatkArgs[ind++] = breakpoints.get(i); } gatkArgs[ind++] =
   * "-o"; gatkArgs[ind++] = tmpOutput + ".bam"; gatkArgs[ind++] = "-l";
   * gatkArgs[ind++] = "ERROR"; if (debug) { gatkArgs[gatkArgs.length - 1] =
   * "INFO"; } RunGatkCommand.run(gatkArgs, debug); FileOutputStream output;
   * output = new FileOutputStream(tmpOutput + ".sam"); RunSystemCommand.run(new
   * ProcessBuilder("samtools", "view", "-h", tmpOutput + ".bam"), output, null,
   * System.out, "samtools", debug); RunSystemCommand.run(new
   * ProcessBuilder("samtools", "view", "-f", "4", bamFile), output, null,
   * System.out, "samtools", debug); output.close(); } else { FileOutputStream
   * output; output = new FileOutputStream(tmpOutput + ".sam");
   * RunSystemCommand.run(new ProcessBuilder("samtools", "view", "-h", "-f",
   * "4", bamFile), output, null, System.out, "samtools", debug);
   * output.close(); } Scanner sc = new Scanner(new File(tmpOutput + ".sam"));
   * PrintWriter writer = new PrintWriter(fqOutputFile); while
   * (sc.hasNextLine()) { String line = sc.nextLine(); if (line.charAt(0) !=
   * '@') { String[] tokens = line.split("\t"); if
   * (BigInteger.valueOf(Long.valueOf(tokens[1]).byteValue()).testBit(6)) {
   * writer.println("@" + tokens[0] + "-1"); } else { writer.println("@" +
   * tokens[0] + "-2"); } if
   * (BigInteger.valueOf(Long.valueOf(tokens[1]).byteValue()).testBit(4)) {
   * writer.println(Genome.rc(tokens[9])); } else { writer.println(tokens[9]); }
   * writer.println("+"); if
   * (BigInteger.valueOf(Long.valueOf(tokens[1]).byteValue()).testBit(4)) {
   * writer.println(Genome.reverse(tokens[10])); } else {
   * writer.println(tokens[10]); } } } sc.close(); writer.close(); new
   * File(tmpOutput + ".bam").delete(); new File(tmpOutput + ".bai").delete();
   * new File(tmpOutput + ".sam").delete(); }
   */

  /*
   * public void localAlign(String fqFile, String faRefFile, String
   * bwaIndexPrefix, String bamOutputFile) throws FileNotFoundException,
   * InterruptedException, IOException, ExecutionException, Exception { String
   * tmpOutput = prefix + ".tmp" + System.currentTimeMillis(); FileOutputStream
   * output = new FileOutputStream(tmpOutput + "-1.sam"); if (debug) {
   * RunSystemCommand.run(new ProcessBuilder("bwa", "bwasw", "-M",
   * bwaIndexPrefix, fqFile), output, null, System.out, "BWA-SW", debug); } else
   * { RunSystemCommand.run(new ProcessBuilder("bwa", "bwasw", "-M",
   * bwaIndexPrefix, fqFile), output, null, null, null, debug); }
   * output.close(); // Sorting must be done before calmd, or else calmd will be
   * very slow addReadGroups(tmpOutput + "-1.sam", tmpOutput + "-2.bam"); output
   * = new FileOutputStream(tmpOutput + "-3.bam"); // calmd is getting stuck! I
   * think the error stream is getting filled up // before it has a chance to be
   * read from. Need concurrent threads? RunSystemCommand.run(new
   * ProcessBuilder("samtools", "calmd", "-b", tmpOutput + "-2.bam", faRefFile),
   * output, null, System.out, "samtools", debug); output.close(); output = new
   * FileOutputStream(bamOutputFile); RunSystemCommand.run(new
   * ProcessBuilder("samtools", "view", "-b", "-F", "256", tmpOutput +
   * "-3.bam"), output, null, System.out, "samtools", debug); output.close();
   * RunSystemCommand.run(new ProcessBuilder("samtools", "index",
   * bamOutputFile), System.out, "samtools", System.out, "samtools", debug); new
   * File(tmpOutput + "-1.sam").delete(); new File(tmpOutput +
   * "-2.bam").delete(); new File(tmpOutput + "-3.bam").delete(); }
   */

  // Helper class to be used by the extractAndRealign method
  public class Region {

    public String chr;
    public long   start;
    public long   end;

    public Region(String chr, long start, long end) {
      this.chr = chr;
      this.start = start;
      this.end = end;
    }

  }

  public void writeSamRecordToFastq(String[] tokens, PrintWriter fqOut) {
    if (BigInteger.valueOf(Long.valueOf(tokens[1]).byteValue()).testBit(6)) {
      fqOut.println("@" + tokens[0] + "-1");
    } else {
      fqOut.println("@" + tokens[0] + "-2");
    }
    if (BigInteger.valueOf(Long.valueOf(tokens[1]).byteValue()).testBit(4)) {
      fqOut.println(Genome.rc(tokens[9]));
      fqOut.println("+");
      fqOut.println(Genome.reverse(tokens[10]));
    } else {
      fqOut.println(tokens[9]);
      fqOut.println("+");
      fqOut.println(tokens[10]);
    }
  }

  public long getAlignmentStart(long pos, String cigar) {
    if (!cigar.matches("([0-9]+(M|I|D|N|S|H|P|=|X))+")) {
      throw new ShearRuntimeException("Invalid CIGAR string encountered: " + cigar);
    }
    long alignmentStart = pos;
    String currentNumStr = "";
    CIGAR_SCAN: for (int i = 0; i < cigar.length(); i++) {
      char c = cigar.charAt(i);
      if (c >= '0' && c <= '9') {
        currentNumStr += c;
      } else {
        long currentNum = Long.valueOf(currentNumStr);
        currentNumStr = "";
        switch (c) {
          case 'H':
            break;
          case 'S':
            alignmentStart -= currentNum;
            break CIGAR_SCAN;
          case 'M':
          case 'I':
          case 'D':
          case 'N':
          case 'P':
          case '=':
          case 'X':
            break CIGAR_SCAN;
          default:
            throw new ShearRuntimeException("Invalid CIGAR string encountered: " + cigar);
        }
      }
    }
    return alignmentStart;
  }

  public long getAlignmentEnd(long pos, String cigar) {
    if (!cigar.matches("([0-9]+(M|I|D|N|S|H|P|=|X))+")) {
      throw new ShearRuntimeException("Invalid CIGAR string encountered: " + cigar);
    }
    long alignmentEnd = pos;
    String currentNumStr = "";
    for (int i = 0; i < cigar.length(); i++) {
      char c = cigar.charAt(i);
      if (c >= '0' && c <= '9') {
        currentNumStr += c;
      } else {
        long currentNum = Long.valueOf(currentNumStr);
        currentNumStr = "";
        switch (c) {
          case 'H':
          case 'I':
          case 'P':
            break;
          case 'S':
            if (pos != alignmentEnd) {
              alignmentEnd += currentNum;
            }
            break;
          case 'M':
          case 'D':
          case 'N':
          case '=':
          case 'X':
            alignmentEnd += currentNum;
            break;
          default:
            throw new ShearRuntimeException("Invalid CIGAR string encountered: " + cigar);
        }
      }
    }
    return alignmentEnd - 1l;
  }

  public void extractAndRealign(
      ArrayList<String> breakpoints, String bamFile, String faRefFile, String bwaIndexPrefix, String bamOutputFile)
      throws FileNotFoundException, IOException, Exception {

    // Parameters
    int extractionPadding = 10;
    int regionPadding = 25;
    int newRegionSeparation = 5;

    // Temp output prefix
    String tmpOutput = prefix + ".tmp" + System.currentTimeMillis();

    // Extract breakpoint and unmapped reads to temp BAM file
    if (breakpoints.size() > 0) {
      String[] gatkArgs = new String[12 + 2 * breakpoints.size()];
      int ind = 0;
      gatkArgs[ind++] = "-T";
      gatkArgs[ind++] = "PrintReads";
      gatkArgs[ind++] = "-I";
      gatkArgs[ind++] = bamFile;
      gatkArgs[ind++] = "-R";
      gatkArgs[ind++] = faRefFile;
      gatkArgs[ind++] = "--interval_padding";
      gatkArgs[ind++] = String.valueOf(extractionPadding);
      for (int i = 0; i < breakpoints.size(); i++) {
        gatkArgs[ind++] = "--intervals";
        gatkArgs[ind++] = breakpoints.get(i);
      }
      gatkArgs[ind++] = "-o";
      gatkArgs[ind++] = tmpOutput + ".extracted.bam";
      gatkArgs[ind++] = "-l";
      gatkArgs[ind++] = "ERROR";
      if (debug) {
        gatkArgs[gatkArgs.length - 1] = "INFO";
      }
      RunGatkCommand.run(gatkArgs, debug);
      FileOutputStream output = new FileOutputStream(tmpOutput + ".extracted.sam");
      RunSystemCommand.run(
          new ProcessBuilder("samtools", "view", "-h", tmpOutput + ".extracted.bam"),
          output,
          null,
          System.out,
          "samtools",
          debug);
      RunSystemCommand.run(
          new ProcessBuilder("samtools", "view", "-f", "4", bamFile),
          output,
          null,
          System.out,
          "samtools",
          debug);
      output.close();

      // Only extract unmapped reads to temp BAM file (if no breakpoints)
    } else {
      FileOutputStream output = new FileOutputStream(tmpOutput + ".extracted.sam");
      RunSystemCommand.run(
          new ProcessBuilder("samtools", "view", "-h", "-f", "4", bamFile),
          output,
          null,
          System.out,
          "samtools",
          debug);
      output.close();
    }

    // Set up scanner and writers
    Scanner sc = new Scanner(new File(tmpOutput + ".extracted.sam"));
    ArrayList<Region> bpRegions = new ArrayList<Region>();
    int regionCtr = 0;
    String currentRegionChr = "";
    long currentRegionStart = 0l;
    long currentRegionEnd = 0l;
    PrintWriter breakpointReadsOut = new PrintWriter(tmpOutput + ".region" + regionCtr + ".fq");
    PrintWriter unmappedReadsWriter = new PrintWriter(tmpOutput + ".unmapped.fq");

    // Loop through temp BAM file and sort reads into FASTQ files
    // Sorts by region reads and unmapped reads
    // Keep track of total number of reads present
    int extractedReads = 0;
    while (sc.hasNextLine()) {
      String line = sc.nextLine();
      if (line.charAt(0) != '@') {
        extractedReads++;
        String[] tokens = line.split("\t");
        if (BigInteger.valueOf(Long.valueOf(tokens[1]).byteValue()).testBit(2)) {
          writeSamRecordToFastq(tokens, unmappedReadsWriter);
        } else {
          String alignmentChr = tokens[2];
          long alignmentStart = getAlignmentStart(Long.valueOf(tokens[3]), tokens[5]);
          long alignmentEnd = getAlignmentEnd(Long.valueOf(tokens[3]), tokens[5]);
          if (!alignmentChr.equals(currentRegionChr) || alignmentStart > currentRegionEnd + newRegionSeparation) {
            bpRegions.add(new Region(currentRegionChr, currentRegionStart, currentRegionEnd));
            breakpointReadsOut.close();
            currentRegionChr = alignmentChr;
            currentRegionStart = alignmentStart;
            currentRegionEnd = alignmentEnd;
            regionCtr++;
            breakpointReadsOut = new PrintWriter(tmpOutput + ".region" + regionCtr + ".fq");
          } else {
            if (alignmentStart < currentRegionStart) {
              currentRegionStart = alignmentStart;
            }
            if (alignmentEnd > currentRegionEnd) {
              currentRegionEnd = alignmentEnd;
            }
          }
          writeSamRecordToFastq(tokens, breakpointReadsOut);
        }
      }
    }

    // Write last region
    bpRegions.add(new Region(currentRegionChr, currentRegionStart, currentRegionEnd));

    // Remove dummy region0 reads
    new File(tmpOutput + ".region0.fq").delete();

    // Close scanner and writers
    sc.close();
    breakpointReadsOut.close();
    unmappedReadsWriter.close();

    // Do local realignment for unmapped reads
    FileOutputStream newAlignmentOut = new FileOutputStream(tmpOutput + ".extracted.merged.sam");
    if (debug) {
      RunSystemCommand.run(
          new ProcessBuilder("bwa", "bwasw", "-M", bwaIndexPrefix, tmpOutput + ".unmapped.fq"),
          newAlignmentOut,
          null,
          System.out,
          "BWA-SW",
          debug);
    } else {
      RunSystemCommand.run(
          new ProcessBuilder("bwa", "bwasw", "-M", bwaIndexPrefix, tmpOutput + ".unmapped.fq"),
          newAlignmentOut,
          null,
          null,
          null,
          debug);
    }

    PrintWriter newAlignmentWriter = new PrintWriter(newAlignmentOut);

    // Do local realignment for breakpoint regions
    Genome gen = new Genome(faRefFile);
    for (int i = 1; i < bpRegions.size(); i++) {

      // Write region
      PrintWriter faOut = new PrintWriter(tmpOutput + ".region" + i + ".fa");
      Region currentRegion = bpRegions.get(i);
      faOut.println(">" + currentRegion.chr);
      long paddedRegionStart = Math.max(1l, currentRegion.start - regionPadding);
      long paddedRegionEnd = Math.min(gen.getRegionLength(currentRegion.chr), currentRegion.end + regionPadding);
      String regionSeq = gen.getSequence(currentRegion.chr, paddedRegionStart, paddedRegionEnd);
      faOut.println(regionSeq);
      faOut.close();

      // Generate BWA index for region
      if (debug) {
        RunSystemCommand.run(
            new ProcessBuilder("bwa", "index", tmpOutput + ".region" + i + ".fa"),
            System.out,
            "BWA",
            System.out,
            "BWA",
            debug);
      } else {
        RunSystemCommand.run(
            new ProcessBuilder("bwa", "index", tmpOutput + ".region" + i + ".fa"),
            null,
            null,
            null,
            null,
            debug);
      }

      // Do local realignment
      FileOutputStream samOut = new FileOutputStream(tmpOutput + ".region" + i + ".sam");
      if (debug) {
        RunSystemCommand.run(new ProcessBuilder("bwa", "bwasw", "-T", "10", "-c", "2", "-M", tmpOutput + ".region" + i
            + ".fa", tmpOutput + ".region" + i + ".fq"), samOut, null, System.out, "BWA-SW", debug);
      } else {
        RunSystemCommand.run(new ProcessBuilder("bwa", "bwasw", "-T", "10", "-c", "2", "-M", tmpOutput + ".region" + i
            + ".fa", tmpOutput + ".region" + i + ".fq"), samOut, null, null, null, debug);
      }
      samOut.close();

      // Merge with main realignment
      sc = new Scanner(new File(tmpOutput + ".region" + i + ".sam"));
      while (sc.hasNextLine()) {
        String line = sc.nextLine();
        if (line.charAt(0) != '@') {
          String[] tokens = line.split("\t");
          for (int k = 0; k < tokens.length; k++) {
            if (k == 3) {
              long adjustedPos = Long.valueOf(tokens[k]) == 0l ? 0l : Long.valueOf(tokens[k]) + paddedRegionStart - 1l;
              newAlignmentWriter.print(adjustedPos + "\t");
            } else if (k == tokens.length - 1) {
              newAlignmentWriter.print(tokens[k] + "\n");
            } else {
              newAlignmentWriter.print(tokens[k] + "\t");
            }
          }
        }
      }
      sc.close();

      // Cleanup temp files
      new File(tmpOutput + ".region" + i + ".fa").delete();
      new File(tmpOutput + ".region" + i + ".fa.amb").delete();
      new File(tmpOutput + ".region" + i + ".fa.ann").delete();
      new File(tmpOutput + ".region" + i + ".fa.bwt").delete();
      new File(tmpOutput + ".region" + i + ".fa.pac").delete();
      new File(tmpOutput + ".region" + i + ".fa.sa").delete();
      new File(tmpOutput + ".region" + i + ".fq").delete();
      new File(tmpOutput + ".region" + i + ".sam").delete();

    }
    gen.close();

    // Close output
    newAlignmentWriter.close();
    newAlignmentOut.close();

    // Process in the usual way
    if (extractedReads > 0) {

      // Remove secondary alignments and unmapped
      // Must remove reads that are secondary AND unmapped alignments
      // before running Picard (i.e. addReadGroupsAndSort)
      newAlignmentOut = new FileOutputStream(tmpOutput + ".extracted.merged.primary-only.bam");
      RunSystemCommand.run(new ProcessBuilder("samtools", "view", "-S", "-b", "-F", "260", tmpOutput
          + ".extracted.merged.sam"), newAlignmentOut, null, System.out, "samtools", debug);
      newAlignmentOut.close();

      // Add read groups and sort
      addReadGroupsAndSort(tmpOutput + ".extracted.merged.primary-only.bam", tmpOutput
          + ".extracted.merged.primary-only.sorted.rg-added.bam");

      // Run calmd on alignment
      newAlignmentOut = new FileOutputStream(bamOutputFile);
      RunSystemCommand.run(
          new ProcessBuilder("samtools", "calmd", "-b", tmpOutput
              + ".extracted.merged.primary-only.sorted.rg-added.bam", faRefFile),
          newAlignmentOut,
          null,
          System.out,
          "samtools",
          debug);
      newAlignmentOut.close();

      // No reads, just output the headers
    } else {

      newAlignmentOut = new FileOutputStream(bamOutputFile);
      RunSystemCommand.run(
          new ProcessBuilder("samtools", "view", "-b", "-S", "-H", tmpOutput + ".extracted.merged.sam"),
          newAlignmentOut,
          null,
          System.out,
          "samtools",
          debug);
      newAlignmentOut.close();

    }

    // Index alignment
    RunSystemCommand.run(
        new ProcessBuilder("samtools", "index", bamOutputFile),
        System.out,
        "samtools",
        System.out,
        "samtools",
        debug);

    // Cleanup temp files
    new File(tmpOutput + ".extracted.bam").delete();
    new File(tmpOutput + ".extracted.bai").delete();
    new File(tmpOutput + ".extracted.sam").delete();
    new File(tmpOutput + ".extracted.merged.sam").delete();
    new File(tmpOutput + ".extracted.merged.primary-only.bam").delete();
    new File(tmpOutput + ".extracted.merged.primary-only.sorted.rg-added.bam").delete();
    new File(tmpOutput + ".unmapped.fq").delete();

  }

  /*
   * // This should be cleaned up! // Will result in duplicate reads (albeit
   * with different names) // When unmapped reads are "re-added" after BWA-SW
   * public void mergeAlignments(String newAlignment, String firstAlignment,
   * String secondAlignment) throws ExecutionException, IOException,
   * InterruptedException { String tmpOutput = prefix + ".tmp" +
   * System.currentTimeMillis(); RunSystemCommand.run(new
   * ProcessBuilder("samtools", "sort", secondAlignment, tmpOutput), System.out,
   * "samtools", System.out, "samtools", debug); RunSystemCommand.run(new
   * ProcessBuilder("samtools", "merge", newAlignment, firstAlignment, tmpOutput
   * + ".bam"), System.out, "samtools", System.out, "samtools", debug);
   * RunSystemCommand.run(new ProcessBuilder("samtools", "index", newAlignment),
   * System.out, "samtools", System.out, "samtools", debug); new File(tmpOutput
   * + ".bam").delete(); }
   */

  public void estimateHeterogeneityAndOutput(
      String predFile, String alignmentFile, String faRefFile, TreeSet<Variant> variants) throws Exception {

    // Estimate heterogeneity for all variants
    System.out.println("[SHEAR] Estimating Heterogeneity...");
    ArrayList<CrestPrediction> crestPreds = getCrestPreds(predFile);
    crestPreds = processCrestPreds(alignmentFile, faRefFile, crestPreds);

    // Filter CREST results
    System.out.println("[SHEAR] Filtering Results...");
    crestPreds = filterCrestPreds(crestPreds);

    // Convert crestPreds to variants, fix this later to just use variants
    for (CrestPrediction pred : crestPreds) {
      if (pred.type.equals("DEL")) {
        Deletion del = new Deletion(pred.bp1region, pred.bp1 + 1, pred.bp2 - pred.bp1 - 1);
        del.setHeterogeneity(pred.heterogeneityPercent);
        variants.add(del);
      } else if (pred.type.equals("INS")) {
        TandemDuplication tdup = new TandemDuplication(pred.bp1region, pred.bp2, pred.bp1 - pred.bp2 + 1);
        tdup.setHeterogeneity(pred.heterogeneityPercent);
        variants.add(tdup);
      } else if (pred.type.equals("CTX") || pred.type.equals("ITX")) {
        Translocation xloc = new Translocation(pred.bp1region, pred.bp1, pred.bp2region, pred.bp2);
        xloc.setHeterogeneity(pred.heterogeneityPercent);
        variants.add(xloc);
        System.out.printf("[SHEAR]       Skipping ITX/CTX for SDI output "
            + "(Just SNPs, INDELs, deletions, inversions, and tandem duplications for now)\n");
      } else if (pred.type.equals("INV") && pred.bp1 < pred.bp2) {
        Inversion inv = new Inversion(pred.bp1region, pred.bp1 + 1, pred.bp2 - pred.bp1 - 1);
        inv.setHeterogeneity(pred.heterogeneityPercent);
        variants.add(inv);
      }
    }

    System.out.println("[SHEAR] Filtering SNPs/INDELs overlapping with SVs...");
    TreeSet<Variant> filteredVariants = new TreeSet<Variant>();
    for (Variant v : variants) {
      if (v instanceof SNP || v.isIndel()) {
        boolean overlaps = false;
        for (Variant sv : variants) {
          if (!(sv instanceof SNP) && !sv.isIndel() && v.adjacentTo(sv)) {
            overlaps = true;
            break;
          }
        }
        if (!overlaps) {
          filteredVariants.add(v);
        }
      } else {
        filteredVariants.add(v);
      }
    }

    if (minHetThresh > 0) {
      System.out.printf("[SHEAR] Filtering variants below %.2f%% estimated heterogeneity level...\n", minHetThresh);
      TreeSet<Variant> unfilteredVariants = filteredVariants;
      filteredVariants = new TreeSet<Variant>();
      for (Variant v : unfilteredVariants) {
        if (v.getHeterogeneity() >= minHetThresh) {
          filteredVariants.add(v);
        }
      }
    }

    System.out.println("[SHEAR] Outputing Results...");
    outputResults(filteredVariants);

    System.out.println("[SHEAR] All Finished!");
  }

  /**
   * Main method to run the ShearSV module.
   * 
   * @throws InterruptedException
   * @throws IOException
   * @throws Exception
   */
  public void execute() throws InterruptedException, IOException, Exception {

    if (hetOnly) {

      estimateHeterogeneityAndOutput(origPreds, origBamFile, origFaRefFile, new TreeSet<Variant>());

    } else {

      // Initialize GFServer
      GFServer gfServer = new GFServer(prefix, origTwoBitRefFile, debug);

      try {

        // Start GFServer
        if (debug) {
          System.out.println("[SHEAR] Using port " + gfServer.getPortNum() + " for gfServer");
          System.out.println("[SHEAR] Starting gfServer...");
        }
        gfServer.start();
        if (debug) {
          System.out.println("[SHEAR] Waiting for gfServer to be ready...");
        }
        boolean gfServerStarted = gfServer.waitForReady();
        while (!gfServerStarted) {
          if (!gfServer.isValidTwoBitFile()) {
            throw new ShearRuntimeException("Not a valid two-bit reference file: " + origTwoBitRefFile);
          }
          if (debug) {
            System.out.println("[SHEAR] Port already in use. Restarting gfServer...");
          }
          gfServer.clearLogFiles();
          gfServer.close();
          gfServer = new GFServer(prefix, origTwoBitRefFile, debug, gfServer.getPortNum());
          if (debug) {
            System.out.println("[SHEAR] Using port " + gfServer.getPortNum() + " for gfServer");
            System.out.println("[SHEAR] Starting gfServer...");
          }
          gfServer.start();
          if (debug) {
            System.out.println("[SHEAR] Waiting for gfServer to be ready...");
          }
          gfServerStarted = gfServer.waitForReady();
        }
        if (debug) {
          System.out.println("[SHEAR] gfServer is ready!");
        }

        // Extract region + unmapped reads, if specified
        // Also extra original read group information
        if (region != null) {
          System.out.println("[SHEAR] Extracting reads from region " + region + "...");
          String tmpOutput = prefix + ".tmp" + System.currentTimeMillis();
          FileOutputStream output = new FileOutputStream(tmpOutput + ".bam");
          if (debug) {
            RunSystemCommand.run(
                new ProcessBuilder("samtools", "view", "-b", origBamFile, region),
                output,
                null,
                System.out,
                "samtools",
                debug);
          } else {
            RunSystemCommand.run(
                new ProcessBuilder("samtools", "view", "-b", origBamFile, region),
                output,
                null,
                null,
                null,
                debug);
          }
          output.close();
          System.out.println("[SHEAR] Extracting original read group information...");
          prepareOrigBamFile(tmpOutput + ".bam");
          new File(tmpOutput + ".bam").delete();
        } else {
          System.out.println("[SHEAR] Extracting original read group information...");
          prepareOrigBamFile(origBamFile);
        }

        // Extract unmapped reads
        // System.out.println("[SHEAR] Extracting Unmapped Reads...");
        // extractRelevantReads(new ArrayList<String>(), prefix +
        // "-orig-with-rg.bam", origFaRefFile, prefix + "." + 0 +
        // ".extracted.fq");

        // Create new alignment with locally-remapped originally-unmapped reads
        // System.out.println("[SHEAR] Realigning using BWA-SW local alignment...");
        // localAlign(prefix + "." + 0 + ".extracted.fq", origFaRefFile,
        // origBwaIndexPrefix, prefix + "." + 0 + ".extracted.aligned.bam");
        // System.out.println("[SHEAR] Merging BWA-SW-aligned reads back to original alignment...");
        // mergeAlignments(prefix + "." + 0 + ".extracted.aligned.merged.bam",
        // prefix + "-orig-with-rg.bam", prefix + "." + 0 +
        // ".extracted.aligned.bam");

        // Initialize set of variants
        TreeSet<Variant> variants = new TreeSet<Variant>();

        // Predict SNPs/INDELs
        if (!svOnly) {
          System.out.println("[SHEAR] Predicting SNPs and INDELs...");
          variants = runHaplotypeCaller(prefix + "-orig-with-rg.bam", origFaRefFile, prefix + ".snpindel.vcf");
        }

        // Initial SV predictions
        String origCrestPredsFile;
        if (origPreds == null) {
          System.out.println("[SHEAR] Making SV Predictions with CREST...");
          runCrest(prefix + "-orig-with-rg.bam", prefix + ".0", origFaRefFile, gfServer, false);
          origCrestPredsFile = prefix + ".0.predSV.txt";
        } else {
          origCrestPredsFile = origPreds;
        }

        // runCrest(prefix + "." + 0 + ".extracted.aligned.merged.bam", prefix +
        // ".0", origFaRefFile, gfServer);

        // Look at breakpoints from SV predictions
        System.out.println("[SHEAR] Extracting Breakpoints...");
        ArrayList<String> breakpoints = getBreakpoints(origCrestPredsFile);
        System.out.println("[SHEAR] Breakpoints:");
        boolean hasNewBreakpoints = true;
        if (breakpoints.size() > 0) {
          for (int i = 0; i < breakpoints.size(); i++) {
            System.out.println("[SHEAR]   " + breakpoints.get(i));
          }
        } else {
          System.out.println("[SHEAR]   None");
          hasNewBreakpoints = false;
        }

        // Continue iterations until no new breakpoints from SV predictions
        // This needs to go at least one iteration in order to generate
        // *.1.predSV.txt, which later code assumes
        int ctr = 0;
        while (ctr < 5 && (ctr == 0 || hasNewBreakpoints)) {
          ctr++;

          // Extract relevant reads and realign
          System.out.println("[SHEAR] Extracting and Realigning Breakpoint and Unmapped Reads...");
          extractAndRealign(breakpoints, prefix + "-orig-with-rg.bam", origFaRefFile, origBwaIndexPrefix, prefix + "."
              + ctr + ".extracted.aligned.bam");

          /*
           * System.out.println(
           * "[SHEAR] Extracting Breakpoint and Unmapped Reads..." );
           * extractRelevantReads(breakpoints, prefix + "-orig-with-rg.bam",
           * origFaRefFile, prefix + "." + ctr + ".extracted.fq");
           * System.out.println
           * ("[SHEAR] Realigning using BWA-SW local alignment...");
           * localAlign(prefix + "." + ctr + ".extracted.fq", origFaRefFile,
           * origBwaIndexPrefix, prefix + "." + ctr + ".extracted.aligned.bam");
           * if (!debug) { new File(prefix + "." + ctr +
           * ".extracted.fq").delete(); }
           */

          // Make new SV predictions
          System.out.println("[SHEAR] Making SV Predictions with CREST...");
          runCrest(
              prefix + "." + ctr + ".extracted.aligned.bam",
              prefix + "." + ctr + "",
              origFaRefFile,
              gfServer,
              true);

          // Look at breakpoints from SV predictions
          System.out.println("[SHEAR] Extracting Breakpoints...");
          ArrayList<String> newBreakpoints = getBreakpoints(prefix + "." + ctr + ".predSV.txt");
          hasNewBreakpoints = false;
          for (int i = 0; i < newBreakpoints.size(); i++) {
            if (!breakpoints.contains(newBreakpoints.get(i))) {
              System.out.println("[SHEAR] New breakpoint: " + newBreakpoints.get(i));
              breakpoints.add(newBreakpoints.get(i));
              hasNewBreakpoints = true;
            }
          }

        }

        // No new breakpoints from SV predictions
        System.out.println("[SHEAR] No new breakpoints detected");
        if (!debug) {
          new File(prefix + "-orig-with-rg.bam").delete();
          new File(prefix + "-orig-with-rg.bam.bai").delete();
        }

        estimateHeterogeneityAndOutput(prefix + "." + ctr + ".predSV.txt", prefix + "." + ctr
            + ".extracted.aligned.bam", origFaRefFile, variants);

        // Cleanup
        if (!debug) {
          for (int i = 0; i <= ctr; i++) {
            new File(prefix + "." + i + ".predSV.txt").delete();
            new File(prefix + "." + i + ".extracted.aligned.bam").delete();
            new File(prefix + "." + i + ".extracted.aligned.bam.bai").delete();
          }
          new File(prefix + ".snpindel.vcf").delete();
          new File(prefix + ".snpindel.vcf.idx").delete();
        }

      } catch (Exception e) {
        throw e;

      } finally {
        if (debug) {
          System.out.println("[SHEAR] Closing gfServer");
        }
        gfServer.close();
        if (!debug) {
          new File(prefix + ".gfServer.log").delete();
          new File(prefix + ".gfServer.out").delete();
        }
      }

    }

  }
}
