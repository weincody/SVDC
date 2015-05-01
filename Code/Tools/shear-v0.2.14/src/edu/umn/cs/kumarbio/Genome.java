/*
 * Sean Landman
 * Data Mining for Biomedical Informatics Group
 * University of Minnesota - Twin Cities
 */

package edu.umn.cs.kumarbio;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Represents a genome, parses FASTA and FAI files, and allows for fast access
 * of sequences.
 */
public class Genome {

  // Access to FASTA file
  private RandomAccessFile         genomeFile;

  // Information from FAI file
  private ArrayList<String>        regionNames;
  private HashMap<String, Long>    regionLengths;
  private HashMap<String, Long>    regionOffsets;
  private HashMap<String, Integer> regionLineLengths;
  private HashMap<String, Integer> regionLineBytes;

  /**
   * Constructs a new Genome object using the path to the FASTA and FAI files
   * (i.e. the part preceding the .fa and .fa.fai extensions). Will parse the
   * FAI file for correctness and then use it to index into the FASTA file for
   * future operations. Throws a FileFoundFoundException for bad paths or a
   * FileFormatException for unexpected formats.
   * 
   * @param path
   *        Path to the FASTA and FAI files, not including file extension.
   * @throws FileNotFoundException
   *         If bad path or missing file(s).
   * @throws FileFormatException
   *         If bad format in the FASTA or FAI files.
   */
  public Genome(String path) throws FileNotFoundException, FileFormatException {
    openFastaFile(path);
    openFastaIndexFile(path + ".fai");
  }

  // Open FASTA file for reading.
  private void openFastaFile(String filePath) throws FileNotFoundException {
    genomeFile = new RandomAccessFile(filePath, "r");
  }

  // Scan in FAI file, checking for format issues.
  private void openFastaIndexFile(String filePath) throws FileNotFoundException, FileFormatException {
    regionNames = new ArrayList<String>();
    regionLengths = new HashMap<String, Long>();
    regionOffsets = new HashMap<String, Long>();
    regionLineLengths = new HashMap<String, Integer>();
    regionLineBytes = new HashMap<String, Integer>();
    int lineCtr = 0;
    try {
      Scanner sc = new Scanner(new File(filePath));
      Pattern p = Pattern.compile("^(\\S+)\\t(\\d+)\\t(\\d+)\\t(\\d+)\\t(\\d+)$");
      while (sc.hasNextLine()) {
        Matcher m = p.matcher(sc.nextLine());
        lineCtr++;
        if (!m.matches()) {
          sc.close();
          throw new FileFormatException("Unexpected format in file <" + filePath + "> at line " + lineCtr);
        }
        String name = m.group(1);
        long length = Long.parseLong(m.group(2));
        long offset = Long.parseLong(m.group(3));
        int lineLength = Integer.parseInt(m.group(4));
        int lineBytes = Integer.parseInt(m.group(5));
        if (lineBytes != (lineLength + 1)) {
          sc.close();
          throw new FileFormatException("Unexpected format in file <" + filePath + "> at line " + lineCtr);
        }
        regionNames.add(name);
        regionLengths.put(name, length);
        regionOffsets.put(name, offset);
        regionLineLengths.put(name, lineLength);
        regionLineBytes.put(name, lineBytes);
      }
      sc.close();
    } catch (NumberFormatException e) {
      throw new FileFormatException("Unexpected format in file <" + filePath + "> at line " + lineCtr, e);
    }
  }

  /**
   * Returns the sequence from this genome in the given region between the
   * specified start and end locations, inclusive. Throws an IOException for
   * problems reading the FASTA file.
   * 
   * @param region
   *        The name of the region containing the desired sequence.
   * @param start
   *        The starting location of the desired sequence (inclusive).
   * @param end
   *        The ending location of the desired sequence (inclusive).
   * @return The sequence string specified by the parameters.
   * @throws IOException
   *         If there are problems reading from the FASTA file.
   */
  public String getSequence(String region, long start, long end) throws IOException {
    if (!regionLengths.containsKey(region)) {
      throw new IllegalArgumentException("Unknown region specified: " + region);
    }
    long length = regionLengths.get(region);
    long offset = regionOffsets.get(region);
    int lineLength = regionLineLengths.get(region);
    int lineBytes = regionLineBytes.get(region);
    if ((end > length) || (start < 1) || (start > end)) {
      throw new IllegalArgumentException("Improper start and/or end locations specified for region " + region + ":"
          + start + "-" + end);
    }
    if (((end - start) + 1) > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Can't retrieve a sequence as big as " + region + ":" + start + "-" + end);
    }
    int desiredSeqLength = (int) ((end - start) + 1);
    StringBuffer seq = new StringBuffer(desiredSeqLength + lineBytes);
    try {
      genomeFile.seek(offset + (((start - 1) / lineLength) * lineBytes) + ((start - 1) % lineLength));
      while (seq.length() != desiredSeqLength) {
        seq.append(genomeFile.readLine());
        if (seq.length() >= desiredSeqLength) {
          return seq.substring(0, desiredSeqLength);
        }
      }
    } catch (IOException e) {
      throw new IOException("Problem reading from the FASTA file", e);
    }
    return null;
  }

  /**
   * Returns a List object containing the region names for this genome.
   * 
   * @return The region names in this genome.
   */
  public List<String> getRegions() {
    return regionNames;
  }

  /**
   * Returns the length of a given region.
   * 
   * @param region
   *        The region to query.
   * @return The length of the query region.
   */
  public long getRegionLength(String region) {
    if (regionLengths.containsKey(region)) {
      return regionLengths.get(region);
    } else {
      throw new IllegalArgumentException("Invalid region specified: " + region);
    }
  }

  /**
   * Closes the genome file by closing the underlying RandomAccessFile object.
   * 
   * @throws IOException
   *         If an I/O error occurs with closing the RandomAccessFile object.
   */
  public void close() throws IOException {
    genomeFile.close();
  }

  /**
   * Returns the reverse of a sequence. Input is a StringBuffer object. Output
   * is a String object.
   * 
   * @param seq
   *        StringBuffer object containing the input string.
   * @return The reverse of the input string.
   */
  public static String reverse(StringBuffer seq) {
    return seq.reverse().toString();
  }

  /**
   * Returns the reverse of a sequence. Input is a String object. Output is a
   * String object.
   * 
   * @param seq
   *        String object containing the input string.
   * @return The reverse of the input string.
   */
  public static String reverse(String seq) {
    return reverse(new StringBuffer(seq));
  }

  /**
   * Returns the reverse complement of a sequence. Input is a StringBuffer
   * object. Output is a String object.
   * 
   * @param seq
   *        StringBuffer object containing the input string.
   * @return The reverse complement of the input string.
   */
  public static String rc(StringBuffer seq) {
    return reverse(complement(seq));
  }

  /**
   * Returns the reverse complement of a sequence. Input is a String object.
   * Output is a String object.
   * 
   * @param seq
   *        String object containing the input string.
   * @return The reverse complement of the input string.
   */
  public static String rc(String seq) {
    return rc(new StringBuffer(seq));
  }

  /**
   * Returns the complement of a sequence. Input is a StringBuffer object.
   * Output is a String object.
   * 
   * @param seq
   *        StringBuffer object containing the input string.
   * @return The complement of the input string.
   */
  public static String complement(StringBuffer seq) {
    StringBuffer complement = new StringBuffer();
    for (int i = 0; i < seq.length(); i++) {
      switch (seq.charAt(i)) {
        case 'A':
          complement.append('T');
          break;
        case 'C':
          complement.append('G');
          break;
        case 'G':
          complement.append('C');
          break;
        case 'T':
          complement.append('A');
          break;
        case 'a':
          complement.append('t');
          break;
        case 'c':
          complement.append('g');
          break;
        case 'g':
          complement.append('c');
          break;
        case 't':
          complement.append('a');
          break;
        default:
          complement.append(seq.charAt(i));
          break;
      }
    }
    return complement.toString();
  }

  /**
   * Returns the complement of a sequence. Input is a String object. Output is a
   * String object.
   * 
   * @param seq
   *        String object containing the input string.
   * @return The complement of the input string.
   */
  public static String complement(String seq) {
    return complement(new StringBuffer(seq));
  }

}
