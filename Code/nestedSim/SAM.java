/*
  Code for comparing corresponding Socrates outputs
*/

import java.util.*;
import java.io.*;
import java.math.*;

//constructor break point
public class SAM {

    public String name;
    public int flags;
    public String ref;
    public int left_aligned;
    public int quality;
    public String cigar;
    public String mates_ref;
    public int mates_left_aligned;
    public int inferred_length;
    public String sequence;
    public String sequence_quality;
    
    public SAM(String line) {
	Scanner scanner = new Scanner(line);
	scanner.useDelimiter("\t");
	try{
	    name = scanner.next();
	    flags = Integer.parseInt(scanner.next());
	    ref = scanner.next();
	    left_aligned = Integer.parseInt(scanner.next());
	    quality = Integer.parseInt(scanner.next());	
	    cigar = scanner.next();
	    mates_ref = scanner.next();
	    mates_left_aligned = scanner.nextInt();
	    inferred_length = scanner.nextInt();
	    sequence = scanner.next();
	    sequence_quality = scanner.next();
	  	    
	} catch (NoSuchElementException e) {
	    System.out.println("Error in parsing line:" + line);
	    System.exit(0);
	}
	scanner.close();
    }

}
