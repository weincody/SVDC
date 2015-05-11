/*
  Code for exporting breakpoint-flanking regions as FASTA file
*/

import java.util.*;
import java.io.*;
import java.math.*;

public class exportAnnotated{
    
    public static String printFastaLine(BP line, int n, char genome) {
	//String fasta = String.format("> %d C1\n%s\n> %d C2\n%s", n, line.C1_realign_consensus, n, line.C2_realign_consensus);
	//can do realign_consensus or anchor_consensus
	String fasta = "";
	if(n != 0){
	    fasta += String.format("\n");
	}
	fasta += String.format("A%d_%s_C1", n, line.C1_anchor);
	if(line.C2_realign_consensus.matches(".*[ATGC].*")) {
	    fasta += String.format("\n");
	    fasta += String.format("A%d_%s_C2", n, line.C2_anchor);
	}
	return fasta;
    }

    public static String printAllFasta(BPs breaks, char genome) {
	String output = "";
	for(int i = 0; i < breaks.size(); i++){
	    output += printFastaLine(breaks.get(i), i, genome);
	}
	return output;
    }
    
    public static void main(String[] args) {
	String fileName = "";
	char genome = '\0';
	if (args.length == 2) {
	    fileName = args[0];
	    genome = args[1].charAt(0);
	} else {
	    System.out.println("more/less than two arguments (file name and S/C)");
	    System.exit(0);
	}
	
	BPs breaks = new BPs();
	breaks.run(fileName);

	String finale = printAllFasta(breaks, genome);
	System.out.println(finale);
	
	//System.out.println(breaks.get(10).C2_realign_chrom);
	//System.out.println(breaks.get(10).C2_realign_pos);
    }
}
