/*
 * Sean Landman Data Mining for Biomedical Informatics Group University of
 * Minnesota - Twin Cities
 */

package edu.umn.cs.kumarbio.shear;

import java.util.Arrays;

/**
 * This is the main SHEAR class. It is used to read the initial command line
 * arguments, launch the appropriate command module, and handle errors.
 */
public class Shear {

  /**
   * Display the usage syntax for SHEAR.
   */
  public static void mainHelp() {
    // Max line length
    // ---------------------------------------------------------------------------------
    System.out.println();
    System.out.println("SHEAR (Sample Heterogeneity Estimation and Assembly by Reference)");
    System.out.println("Version:  0.2.14 (beta)");
    System.out.println("Website:  http://vk.cs.umn.edu/SHEAR/");
    System.out.println();
    System.out.println("Usage:");
    System.out.println("    java -jar SHEAR.jar <command> [options]");
    System.out.println();
    System.out.println("commands:");
    System.out.println("    sv");
    System.out.println("    assemble");
    System.out.println();
    System.out.println("sv options:");
    System.out.println("    -p <prefix>                        Prefix for all generated files. Required.");
    System.out.println();
    System.out.println("    -b <bam_file>                      BAM alignment file containing the input");
    System.out.println("                                       sequences to the assembly. Must have BAM");
    System.out.println("                                       index file (i.e. *.bam.bai) in same");
    System.out.println("                                       directory. Required.");
    System.out.println();
    System.out.println("    -f <fasta_reference_file>          Reference file (FASTA format) to be used");
    System.out.println("                                       for guiding assembly. Must have FASTA");
    System.out.println("                                       index file (i.e. *.fa.fai) and dictionary");
    System.out.println("                                       file (i.e. *.dict) in same directory.");
    System.out.println("                                       Required.");
    System.out.println();
    System.out.println("    -i <bwa_index_prefix>              Prefix of the BWA index files for this");
    System.out.println("                                       reference (i.e. <bwa_index_prefix>.sa,");
    System.out.println("                                       <bwa_index_prefix>.pac, etc.). Uses");
    System.out.println("                                       <fasta_reference_file> by default.");
    System.out.println();
    System.out.println("    -t <two_bit_reference_file>        Reference file (2bit format) to be used");
    System.out.println("                                       for guiding assembly. Required.");
    System.out.println();
    System.out.println("    -r <region>                        Region of the input alignment to analyze");
    System.out.println("                                       (i.e. \"chr1\" for full chromosome or");
    System.out.println("                                       \"chr1:10000-15000\" for a particular");
    System.out.println("                                       region, including endpoints). Coordinates");
    System.out.println("                                       are 1-based. Optional.");
    System.out.println();
    System.out.println("    --preds <crest_pred_file>          Pre-existing *.predSV.txt file generated");
    System.out.println("                                       by CREST for this alignment. Can be");
    System.out.println("                                       useful to use if CREST has already been");
    System.out.println("                                       run on this alignment independently, or");
    System.out.println("                                       if SHEAR was interrupted and you want to");
    System.out.println("                                       continue from the previous intermediary");
    System.out.println("                                       results. Optional.");
    System.out.println();
    System.out.println("    --min-het <min_het_threshold>      Minimum estimated heterogeneity level for");
    System.out.println("                                       a variant to be included in the results.");
    System.out.println("                                       Must be a numeric value between 0 and 1");
    System.out.println("                                       (i.e. \"--min-het 0.2\" would only output");
    System.out.println("                                       variants that have an estimated).");
    System.out.println("                                       heterogeneity level of at least 20%.");
    System.out.println("                                       Optional.");
    System.out.println();
    System.out.println("    --sv-only                          \"SV Only\" prediction mode. Will skip");
    System.out.println("                                       SNP/INDEL prediction and only return");
    System.out.println("                                       *.sdi and *.report results for SVs.");
    System.out.println();
    System.out.println("    -d                                 Debugging mode. Will keep all");
    System.out.println("                                       intermediary files and produce more");
    System.out.println("                                       verbose output information.");
    System.out.println();
    System.out.println("assemble options:");
    System.out.println("    -s <sdi_file>                      SDI file produced by SHEAR's 'sv' command");
    System.out.println("                                       containing the SVs to use to create the");
    System.out.println("                                       new genomic sequence.");
    System.out.println();
    System.out.println("    -f <fasta_reference_file>          Reference file (FASTA format) to be used");
    System.out.println("                                       for guiding assembly.");
    System.out.println();
    System.out.println("    -o <output_fasta_file>             FASTA file to output. Required.");
    System.out.println();
    System.out.println("    -d                                 Debugging mode. Will keep all");
    System.out.println("                                       intermediary files and produce more");
    System.out.println("                                       verbose output information.");
    System.out.println();
    System.exit(0);
    // Max line length
    // ---------------------------------------------------------------------------------
  }

  /**
   * Display a usage error message and exit. Used when the user does not specify
   * proper arguments.
   * 
   * @param message
   *          The message to display to the user regarding what went wrong
   */
  public static void usageError(String message) {
    System.out.println("[SHEAR] !----- ERROR -----!");
    System.out.println("[SHEAR] " + message);
    System.out.println("[SHEAR] SHEAR usage syntax:");
    mainHelp();
    System.exit(1);
  }

  /**
   * Display a runtime error message and exit. Used when there is some
   * unexpected problem with SHEAR.
   * 
   * @param e
   *          Exception that triggered the runtime error
   */
  public static void runtimeError(Exception e) {
    System.out.println("[SHEAR] !----- ERROR -----!");
    System.out.println("[SHEAR] SHEAR has encountered an unexpected problem");
    System.out.println("[SHEAR] Please contact the developers if the problem persists");
    e.printStackTrace();
    System.out.println("[SHEAR] Exiting...");
    System.exit(1);
  }

  /**
   * Main method to check for command, pass the rest of the arguments to the
   * proper command module, and handle error catching.
   * 
   * @param args
   *          Command line arguments passed to SHEAR
   */
  public static void main(String[] args) {

    try {

      // Check for existence of command
      if (args.length == 0) {
        throw new ShearUsageException("Must specify a command");

        // Run SHEAR-SV
      } else if (args[0].equals("sv")) {
        ShearSV shearSV = new ShearSV(Arrays.copyOfRange(args, 1, args.length));
        shearSV.execute();

        // Run SHEAR-Assemble
      } else if (args[0].equals("assemble")) {
        ShearAssemble shearAssemble = new ShearAssemble(Arrays.copyOfRange(args, 1, args.length));
        shearAssemble.execute();

        // Unrecognized command
      } else {
        throw new ShearUsageException("Unrecognized command: " + args[0]);
      }

      // Catch SHEAR runtime exception
    } catch (ShearRuntimeException e) {
      Shear.runtimeError(e);

      // Catch SHEAR usage exception
    } catch (ShearUsageException e) {
      Shear.usageError(e.getMessage());

      // Catch other exceptions
      // Shouldn't happen
    } catch (Exception e) {
      Shear.runtimeError(e);
    }

  }

}
