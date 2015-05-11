/*
  Code for comparing corresponding Socrates outputs
*/

import java.util.*;
import java.io.*;
import java.math.*;

//constructor break point
public class hitTest {

    public char cancer_or_somatic;
    public int ID;
    public int index;
    //public char cancer_or_somatic2;
    //public int index2;

    public hitTest(String line) {
	Scanner scanner = new Scanner(line);
	scanner.useDelimiter("\t");
	try{
	    scanner.useDelimiter("_");
	    String temp = scanner.next();
	    //System.out.println(temp);
	    cancer_or_somatic = temp.charAt(0);
	    ID = Integer.parseInt(temp.substring(1));
	    
	    //Integer.parseInt(temp.replaceAll("[\\D]", ""));
	    //System.out.println(""+ID);
	    scanner.useDelimiter(":");
	    scanner.next();
	    scanner.useDelimiter("_");
	    index = Integer.parseInt(scanner.next().replaceAll("[\\D]", ""));
	    
	    /*
	    scanner.useDelimiter("\t");
	    scanner.next();
	    if scanner.next() != null;
	    scanner.useDelimiter(":");
	     dont have second index for unmatched list
	    try{
	    	cancer_or_somatic1 = scanner.next().charAt(0);
	    } catch (Exception e){
	    	cancer_or_somatic1 = 'N';		
	    }
	    scanner.useDelimiter("_");
	    try{
	    	index1 = Integer.parseInt(s.replaceAll("[\\D]", ""));
	    } catch (Exception e){
	    	index2 = -1;
	    }
	    */
	    
	    /*
	    System.out.println("hitTest cancer_or_somatic: " + cancer_or_somatic);
	    System.out.println("hitTest ID: " + Integer.toString(ID));
	    System.out.println("hitTest index: " + Integer.toString(index));
	    */
	} catch (NoSuchElementException e) {
	    System.out.println("Error in parsing line:" + line);
	    System.exit(0);
	}
	scanner.close();
    }
}
