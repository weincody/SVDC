/*
  Code for comparing corresponding Socrates outputs
*/

import java.util.*;
import java.io.*;
import java.math.*;

public class BPs {
    public ArrayList<BP> book;

    public BPs(){
	book = new ArrayList<BP>();
    }

    public void add(BP break0){
	boolean added = false;
	if(book.size() == 0){
	    book.add(break0);
	    added = true;
	    return;
	}
	book.add(break0);
    }

    public BP get(int i){
	return book.get(i);
    }

    public int size(){
        return book.size();
    }

    public void printBook(){
        for(int i = 0; i < book.size(); i++){
            System.out.println(book.get(i));
        }
    }

    public void run(String f) {
	Scanner scanner = null;

	// attempt to load file
	try {
	    scanner = new Scanner(new File(f), "UTF-8");
	} catch (NullPointerException break0) {
	    System.out.print("Bad file name.");
	    System.exit(0);
	} catch (java.io.FileNotFoundException break0) {
	    System.out.println("File " + f + " not found.");
	    System.exit(0);
	}
	scanner.nextLine();
	while(scanner.hasNextLine()) {
	    BP break0 = new BP(scanner.nextLine());
	    this.add(break0);
	}
    }

    public static void main(String[] args) {
	String fileName = "";
	if (args.length == 1) {
	    fileName = args[0];
	} else {
	    System.out.println("more than one argument");
	    System.exit(0);
	}
	
	BPs breaks = new BPs();
	breaks.run(fileName);
	System.out.println(breaks.get(10).C2_realign_chrom);
	System.out.println(breaks.get(10).C2_realign_pos);
    }
    
}
