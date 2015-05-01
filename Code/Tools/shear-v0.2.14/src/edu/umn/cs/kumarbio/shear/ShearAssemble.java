
package edu.umn.cs.kumarbio.shear;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

import edu.umn.cs.kumarbio.FileFormatException;

public class ShearAssemble {

  private String  sdiFile;
  private String  origFaRefFile;
  private String  newFaRefFile;
  private boolean debug;

  public ShearAssemble(String[] args) throws ShearUsageException {
    debug = false;
    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-s")) {
        sdiFile = args[++i];
      } else if (args[i].equals("-f")) {
        origFaRefFile = args[++i];
      } else if (args[i].equals("-o")) {
        newFaRefFile = args[++i];
      } else if (args[i].equals("-d")) {
        debug = true;
      } else {
        throw new ShearUsageException("Unknown argument: " + args[i]);
      }
    }
    if (sdiFile == null) {
      throw new ShearUsageException("Must provide an SDI file with the '-s' argument");
    }
    if (origFaRefFile == null) {
      throw new ShearUsageException("Must provide a FASTA reference file with the '-f' argument");
    }
    if (newFaRefFile == null) {
      throw new ShearUsageException("Must provide an output FASTA file with the '-o' argument");
    }
    if (!(new File(sdiFile).exists())) {
      throw new ShearUsageException("SDI file does not exist: " + sdiFile);
    }
    if (!(new File(origFaRefFile).exists())) {
      throw new ShearUsageException("FASTA reference file does not exist: " + origFaRefFile);
    }
    System.out.println("[SHEAR] Starting SHEAR-Assemble");
    System.out.println("[SHEAR] SDI file:               " + sdiFile);
    System.out.println("[SHEAR] FASTA reference file:   " + origFaRefFile);
    System.out.println("[SHEAR] FASTA output file:      " + newFaRefFile);
    if (debug) {
      System.out.println("[SHEAR] Debugging mode:         ON");
    } else {
      System.out.println("[SHEAR] Debugging mode:         OFF");
    }
  }

  private SdiLine getNextMutation(Scanner sc, int lineNum) throws FileFormatException {
    if (sc.hasNextLine()) {
      String[] tokens = sc.nextLine().split("\t");
      String region = tokens[0];
      int loc = Integer.valueOf(tokens[1]);
      int diff = Integer.valueOf(tokens[2]);
      String origBases = tokens[3].equals("-") ? "" : tokens[3];
      String newBases = tokens[4].equals("-") ? "" : tokens[4];
      if ((newBases.length() - origBases.length()) != diff) {
        throw new FileFormatException("Unexpected format in file <" + sdiFile + "> at line " + lineNum);
      }
      return new SdiLine(region, loc, origBases, newBases);
    } else {
      return null;
    }
  }

  private boolean mutationPresent(SdiLine nextMutation, String currentRegion, int regionBasesRead) {
    return (nextMutation != null) && nextMutation.region.equals(currentRegion) && (nextMutation.loc <= regionBasesRead);
  }

  public void execute() throws FileFormatException, IOException {

    // Set up I/O
    BufferedReader origFaRefIn = new BufferedReader(new FileReader(origFaRefFile));
    Scanner sdiIn = new Scanner(new File(sdiFile));
    FileWriter newFastaOut = new FileWriter(newFaRefFile);

    // Initialize variables
    String currentRegion = null;
    int regionBasesRead = 0;
    int processedAndAssembledBases = 0;
    StringBuffer toWrite = new StringBuffer(100000);
    int sdiLineCtr = 1;
    SdiLine nextMutation = getNextMutation(sdiIn, sdiLineCtr++);

    // Loop through each line in original FASTA file
    String line = origFaRefIn.readLine();
    while (line != null) {

      // Check for new region
      if (line.charAt(0) == '>') {
        if (toWrite.length() > 0) {
          newFastaOut.write(toWrite + "\n");
          toWrite.delete(0, toWrite.length());
        }
        Scanner sequenceHeaderLineScanner = new Scanner(line);
        currentRegion = sequenceHeaderLineScanner.next().substring(1);
        sequenceHeaderLineScanner.close();
        regionBasesRead = 0;
        processedAndAssembledBases = 0;
        System.out.println("[SHEAR] Parsing region \"" + currentRegion + "\"");
        newFastaOut.write(line + "\n");

        // Process regular (non-header) line
      } else {

        // Add bases to toWrite queue
        toWrite.append(line);
        regionBasesRead += line.length();

        // Check for any relevant mutations for queue
        while (mutationPresent(nextMutation, currentRegion, regionBasesRead)) {
          int diff = nextMutation.newBases.length() - nextMutation.origBases.length();
          System.out.println("[SHEAR]   Processing mutation: " + nextMutation.region + ":" + nextMutation.loc
              + "\tbp change: " + diff);
          if (nextMutation.loc <= processedAndAssembledBases) {
            System.out.println("[SHEAR]     Skipping this mutation. It overlaps with a previously processed mutation.");
          } else {
            while (regionBasesRead < ((nextMutation.loc + nextMutation.origBases.length()) - 1)) {
              line = origFaRefIn.readLine();
              toWrite.append(line);
              regionBasesRead += line.length();
            }
            int i1 = (toWrite.length() - 1) - (regionBasesRead - nextMutation.loc);
            int i2 = ((toWrite.length() - 1) - (regionBasesRead - nextMutation.loc)) + nextMutation.origBases.length();
            if (!toWrite.substring(i1, i2).equals(nextMutation.origBases)) {
              origFaRefIn.close();
              sdiIn.close();
              newFastaOut.close();
              throw new FileFormatException("SDI file does not match reference sequence in file <" + sdiFile
                  + "> at line " + sdiLineCtr++);
            }
            toWrite.replace(i1, i2, nextMutation.newBases);
            processedAndAssembledBases = nextMutation.loc + nextMutation.origBases.length() - 1;
          }
          nextMutation = getNextMutation(sdiIn, sdiLineCtr++);
        }

        // Write any lines to output
        while (toWrite.length() >= 60) {
          newFastaOut.write(toWrite.substring(0, 60) + "\n");
          toWrite.delete(0, 60);
        }

      }

      // Grab the next line
      line = origFaRefIn.readLine();

    }

    // Write anything remaining
    if (toWrite.length() > 0) {
      newFastaOut.write(toWrite + "\n");
      toWrite.delete(0, toWrite.length());
    }

    // Close I/O
    System.out.println("[SHEAR] All Finished!");
    origFaRefIn.close();
    sdiIn.close();
    newFastaOut.close();

  }
}
