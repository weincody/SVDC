/*
  Code for exporting breakpoint-flanking regions as FASTA file
  USAGE: javac hitTest.java hitList.java runResults.java
  java runResults chr11nested > DCresults.txt
  
  OUTPUT:
  insertion
  deletion
  inversion
  false positives
  
*/

import java.util.*;
import java.io.*;
import java.math.*;
import org.apache.commons.lang3.ArrayUtils;

public class runResults{
    

    private static int readInt(BufferedInputStream in) throws IOException {
    int ret = 0;
    boolean dig = false;
    int sf = 0;
    for (int c = 0; (c = in.read()) != -1; ) {
	//System.out.println( "c: " + c);
        if ((c >= '0' && c <= '9')) {
            dig = true;
            ret = ret * 10 + c - '0';
        }
	else if (c == '-'){
	    dig = true;
	    sf = 1;
	}
	else if (dig) break;
    }

    if(dig == false){
	return -2;
    }

    //System.out.println("ret: " + ret);
    if(sf == 1){
	return -1;
    }
    //System.out.println( "ret: " + ret);
    return ret;

    }


    public static int[][] loadInds(String simStem){

	
	BufferedInputStream input1 = null;
	try{
	    input1 = new BufferedInputStream(new FileInputStream( simStem+"_INDS.txt" ) );
	} catch(Exception e){}
	

	int counter = 0;
	try {
	    while (readInt(input1) > -2){
		counter++;
	    }
	} catch ( IOException  eof ) {}

	try {
	    input1.close();
	} catch ( IOException  eof ) {}
	

	int columns = counter/3;
	int rows = 3;
	
	
	//System.out.println( "The length is: " + columns);

	int[][] a = new int[rows][columns];
	BufferedInputStream bis = null;
	try{
	    bis = new BufferedInputStream(new FileInputStream(simStem+"_INDS.txt"));
	} catch(Exception e) {}

	for (int j = 0; j < columns; j++) {
	    for (int i = 0; i < rows; i++) {
		try{
		    a[i][j] = readInt(bis);
		    //System.out.println(a[i][j]);
		} catch(Exception e) {}		
	    }
	    //System.out.println(a[0][j]);
	}

	//System.out.println("runResults.loadInds a[0][columns - 1]: " + Integer.toString(a[0][columns - 1]));

	return a;

    }

    
    public static int[][] loadEvents(String simStem){

	BufferedInputStream input1 = null;
	try{
	    input1 = new BufferedInputStream(new FileInputStream( simStem+"_EVENTS.txt" ) );
	} catch(Exception e){}

	int counter = 0;

	try {
	    while (readInt(input1) > -2){
		counter++;
	    }
	} catch ( IOException  eof ) {}

	try {
	    input1.close();
	} catch ( IOException  eof ) {}

	int columns = counter/31;
	int rows = 31;

	//System.out.println( "The events length is: " + columns);

	int[][] a = new int[rows][columns];
	
	BufferedInputStream bis = null;
	try{
	    bis = new BufferedInputStream(new FileInputStream(simStem+"_EVENTS.txt"));
	} catch(Exception e) {}
	
	for (int i = 0; i < rows; i++) {
	    for (int j = 0; j < columns; j++) {
		try{
		    a[i][j] = readInt(bis);
		} catch(Exception e) {}		
	    }
	}
	
	//System.out.println("runResults.loadInds events[0][columns - 1]: " + Integer.toString(a[0][columns - 1]));

	return a;

    }

    public static void checkHits(hitList matched, int[][] inds, int[][] events){
	
	// insertion, deletion, inversion
	int[] hits = {0, 0, 0};
	// false positive count
	int FP = 0;

	int i;
	int index = 0;
	int[] toCheck = new int[22];
	//System.out.println("toCheck INITIALIZED: " + toCheck[5]);
	int k;
	int j;
	int flag_ev;
	int flag_noev;
	int hit = 0;
	int lastHit = 0;
	
	Integer holder0[] = new Integer[inds[0].length];
	for (index = 0; index < inds[0].length; index++){
	    holder0[index] = (inds[0][index]);
	    //System.out.println("inds[0][index]: " + inds[0][index]);
	    //System.out.println("holder0[index]: " + holder0.get(index));
	}
	Integer holder1[] = new Integer[inds[1].length];
	for (index = 0; index < inds[1].length; index++){
	    holder1[index] = (inds[1][index]);
	}
	Integer holder2[] = new Integer[inds[2].length];
	for (index = 0; index < inds[2].length; index++){
	    holder2[index] = (inds[2][index]);
	}
	
	for(i = 0; i < matched.size(); i++){
	    toCheck = new int[22];
	    j = 0;

	    // tumor indices
	    if(matched.get(i).cancer_or_somatic == 'C'){
		for(j = -10; j < 11; j++){
		    toCheck[j + 10] = Arrays.asList(holder0).indexOf((Integer) (matched.get(i).index + j));
		    //ArrayUtils.find(holder, matched.get(i).index + j);
		    //System.out.println("matched.get(i).index: " + matched.get(i).index);
		    //System.out.println("toCheck C: " + toCheck[j + 10]);
		}
	    }

	    // normal indices
	    else if(matched.get(i).cancer_or_somatic == 'S'){
		for(j = -10; j < 11; j++){
		    toCheck[j + 10] = Arrays.asList(holder1).indexOf(matched.get(i).index + j);
		    //System.out.println("matched.get(i).index: " + matched.get(i).index);
		    //System.out.println("toCheck S: " + toCheck[j + 10]);
		}
	    }
	    // original indices for just Socrates output
	    else{
		for(j = -10; j < 11; j++){
		    toCheck[j + 10] = Arrays.asList(holder2).indexOf(matched.get(i).index + j);
		    //System.out.println("matched.get(i).index: " + matched.get(i).index);
		    //System.out.println("toCheck A: " + toCheck[j + 10]);
		}   
	    }

	    hit = 0;
	    for(k = 0; (k < 31) & (hit == 0); k++){
		hit = 0;
		flag_ev = 0;
		flag_noev = 0;
		for(j = 0; j < 22; j++){
		    //System.out.println(toCheck[j]);
		    if(toCheck[j] != -1){
			//System.out.println("events[k][toCheck[j]]: " + events[k][toCheck[j]]);
			if(events[k][toCheck[j]] > 1){
			    System.out.println("toCheck j: " + toCheck[j]);
			    System.out.println("events[k][toCheck[j]]: " + events[k][toCheck[j]]);
			    flag_ev = events[k][toCheck[j]];
			}
			else{
			    flag_noev = 1;
			}
		    }
		}
		// if hit both event and non, are breakpoint
		if ((flag_ev > 1) & (flag_noev == 1)){
		    hit = flag_ev;
		    break;
		}
		
	    }

	    // case when second map of same breakpoint
	    if((i > 0) && (matched.get(i).ID == matched.get(i - 1).ID) && (matched.get(i).cancer_or_somatic == matched.get(i - 1).cancer_or_somatic)){
		// if already called as a hit, continue
		if(lastHit > 1){}
		else{
		    // if 2nd detected as hit but not first, undo FP and add hit
		    if(hit > 1){
			hits[hit - 2] += 1;
			FP -= 1;			
		    }
		    // if neither were hits, continue
		}
	    }
	    
	    else{
		if(hit > 1){
		    System.out.println("hit: " + hit);
		    hits[hit - 2] += 1;
		}
		else{
		    FP += 1;
		}
	    }
	    lastHit = hit;
	}

	System.out.println(hits[0]);
	System.out.println(hits[1]);
	System.out.println(hits[2]);
	System.out.println(FP);
	System.out.println("\n");
    }
        
    public static void main(String[] args) {
	String simStem = "";
	if (args.length == 1) {
	    simStem = args[0];
	} else {
	    System.out.println("more/less than one argument (simulation root)");
	    System.exit(0);
	}
	
	hitList matched = new hitList();
	matched.run("matchedBPs.txt");

	hitList unmatched = new hitList();
	unmatched.run("unmatchedBPs.txt");
	//System.out.println("annotated");
	hitList annotated = new hitList();
	annotated.run("annotatedBPs.txt");

	int[][] inds = loadInds(simStem);
	int[][] events = loadEvents(simStem);

	//for(int i = 0; i < inds[0].length; i++){
	    //System.out.println(inds[0][i]);
	    //}

	checkHits(matched, inds, events);
	checkHits(unmatched, inds, events);
	checkHits(annotated, inds, events);
	
    }
}
