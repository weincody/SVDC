/*
  Code for comparing corresponding Socrates outputs
*/

import java.util.*;
import java.io.*;
import java.math.*;

//constructor break point
public class BP {

    public int C1_realign_chrom;
    public int C1_realign_pos;
    public char C1_realign_dir;
    public String C1_realign_consensus;
    public String C1_anchor;
    public char C1_anchor_dir;
    public String C1_anchor_consensus;
    public int C1_long_support;
    public int C1_long_support_bases;	
    public int C1_short_support;
    public int C1_short_support_bases;
    public int C1_short_support_max_len;
    public float C1_avg_realign_mapq;

    public int C2_realign_chrom;
    public int C2_realign_pos;
    public char C2_realign_dir;
    public String C2_realign_consensus;
    public String C2_anchor;
    public char C2_anchor_dir;
    public String C2_anchor_consensus;
    public int C2_long_support;
    public int C2_long_support_bases;	
    public int C2_short_support;
    public int C2_short_support_bases;
    public int C2_short_support_max_len;
    public float C2_avg_realign_mapq;
    
    public String BP_condition;

    public BP(String line) {
	Scanner scanner = new Scanner(line);
	scanner.useDelimiter("\t");
	try{
	    scanner.useDelimiter(":");
	    String s = scanner.next();
	    C1_realign_chrom = Integer.parseInt(s.replaceAll("[\\D]", ""));
	    scanner.useDelimiter("\t");
	    C1_realign_pos = Integer.parseInt(scanner.next().substring(1));
	    C1_realign_dir = scanner.next().charAt(0);
	    C1_realign_consensus = scanner.next();
	    C1_anchor = scanner.next();	
	    C1_anchor_dir = scanner.next().charAt(0);
	    C1_anchor_consensus = scanner.next();
	    C1_long_support	= scanner.nextInt();
	    C1_long_support_bases = scanner.nextInt();
	    C1_short_support = scanner.nextInt();
	    C1_short_support_bases = scanner.nextInt();
	    C1_short_support_max_len = scanner.nextInt();
	    C1_avg_realign_mapq = scanner.nextFloat();
	    
	    scanner.useDelimiter(":");
	    String t = scanner.next().replaceAll("[\\D]", "");
	    C2_realign_chrom = Integer.parseInt(t.substring(1));
	    scanner.useDelimiter("\t");
	    String u = scanner.next().replaceAll("[\\D]", "");
	    C2_realign_pos = Integer.parseInt(u.substring(1));
	    C2_realign_dir = scanner.next().charAt(0);
	    C2_realign_consensus = scanner.next();
	    C2_anchor = scanner.next();	
	    C2_anchor_dir = scanner.next().charAt(0);
	    C2_anchor_consensus = scanner.next();
	    C2_long_support	= scanner.nextInt();
	    C2_long_support_bases = scanner.nextInt();
	    C2_short_support = scanner.nextInt();
	    C2_short_support_bases = scanner.nextInt();
	    C2_short_support_max_len = scanner.nextInt();
	    C2_avg_realign_mapq = scanner.nextFloat();	
	    
	    BP_condition = scanner.next();
	    
	} catch (NoSuchElementException e) {
	    System.out.println("Error in parsing line:" + line);
	    System.exit(0);
	}
	scanner.close();
    }
}
