/*
  Code for exporting final matched and unmatched breakpoints
*/

import java.util.*;
import java.io.*;
import java.math.*;

public class exportMatches {
    
    static public String printMatch(SAM book){        
	String bookString = String.format("%s\t%d\t%s\n", book.name, book.flags, book.ref);
	return bookString;
    }

    static public String printMatchOpp(SAM book){        
	String bookString = String.format("%s\t%d\t%s\n", book.ref, book.flags, book.name);
	return bookString;
    }

    static public String printNoMatch(SAM book){
	String bookString = String.format("%s\t%d\n", book.name, book.flags);
	return bookString;
    }

    static public boolean checkMatch(SAMs breaks, int i){
	if(breaks.get(i).ref.matches("\\*")){
	    return false;
	}
	else{
	    return true;
	}
    }

    static public boolean checkMatchOpp(SAMs breaks, SAMs breaksOpp, int i){
	String name = breaks.get(i).name;
	for(int j = 0; j < breaksOpp.size(); j++){
	    if(breaksOpp.get(j).ref == name){
		return true;
	    }
	}
	return false;
    }

    public static ArrayList<String> checkMatches(SAMs breaksSom, SAMs breaksCanc) {
	String matched = "";
	String unmatched = "";
	for(int i = 0; i < breaksSom.size(); i++){
	    if(checkMatch(breaksSom, i)){
		matched += printMatch(breaksSom.get(i));
	    }
	    else if(checkMatchOpp(breaksSom, breaksCanc, i)){
		break;
	    }
	    else{
		unmatched += printNoMatch(breaksSom.get(i));
	    }

	}

	for(int j = 0; j < breaksCanc.size(); j++){
	    if(checkMatch(breaksCanc, j)){
		matched += printMatch(breaksCanc.get(j));
	    }
	    else if(checkMatchOpp(breaksCanc, breaksSom, j)){
		break;
	    }
	    else{
		unmatched += printNoMatch(breaksCanc.get(j));
	    }
	}

	ArrayList<String> output = new ArrayList<String>();
	output.add(matched);
	output.add(unmatched);

	return output;
    }

    public static void exportFiles(String matched, String unmatched){

	try {
	    FileWriter writerM = new FileWriter("matchedBPs.txt");
	    writerM.write(matched);
	    writerM.close();
        } catch(IOException fnfe) { 
            System.out.println(fnfe.getMessage());
        } 

	try {
	    FileWriter writerM = new FileWriter("unmatchedBPs.txt");
	    writerM.write(unmatched);
	    writerM.close();
        } catch(IOException fnfe2) { 
            System.out.println(fnfe2.getMessage());
        } 

    }
    
    public static void main(String[] args) {
	String fileNameSom = "";
	String fileNameCanc = "";
	if (args.length == 2) {
	    fileNameSom = args[0];
	    fileNameCanc = args[1];
	} else {
	    System.out.println("more/less than two arguments (Somatic and Cancer)");
	    System.exit(0);
	}
	
	SAMs breaksSom = new SAMs();
	breaksSom.run(fileNameSom);

	SAMs breaksCanc = new SAMs();
	breaksCanc.run(fileNameCanc);

	ArrayList<String> matchQuery = checkMatches(breaksSom, breaksCanc);
	String matched = matchQuery.get(0);
	String unmatched = matchQuery.get(1);
	
	exportFiles(matched, unmatched);
    }
}
