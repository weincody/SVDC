/*
  Code for comparing corresponding Socrates outputs
*/

import java.util.*;
import java.io.*;
import java.math.*;

public class SAMs {
    public ArrayList<SAM> book;

    public SAMs(){
	book = new ArrayList<SAM>();
    }

    public void add(SAM break0){
	boolean added = false;
	if(book.size() == 0){
	    book.add(break0);
	    added = true;
	    return;
	}
	book.add(break0);
    }

    public SAM get(int i){
	return book.get(i);
    }

    public int size(){
        return book.size();
    }

    public String printSimpleSAMs(){
	String bookString = "";
        for(int i = 0; i < book.size(); i++){
	    bookString += String.format("%s\t%d\t%s\t%s\n", book.get(i).name, book.get(i).flags, book.get(i).ref, book.get(i).left_aligned);
        }
	return bookString;
    }

    public void run(String f) {
	Scanner scanner = null;

	// attempt to load file (count lines)
	try {
	    scanner = new Scanner(new File(f), "UTF-8");
	} catch (NullPointerException break0) {
	    System.out.print("Bad file name.");
	    System.exit(0);
	} catch (java.io.FileNotFoundException break0) {
	    System.out.println("File " + f + " not found.");
	    System.exit(0);
	}
	int count = 0;
	while(scanner.next().charAt(0) == '@'){
	    scanner.nextLine();
	    count++;
	}

	// reload file again (for real)
	scanner = null;
	try {
	    scanner = new Scanner(new File(f), "UTF-8");
	} catch (NullPointerException break0) {
	    System.out.print("Bad file name.");
	    System.exit(0);
	} catch (java.io.FileNotFoundException break0) {
	    System.out.println("File " + f + " not found.");
	    System.exit(0);
	}

	for(int i = 0; i < count; i++) {
	    scanner.nextLine();
	}
	//System.out.println(scanner.nextLine());
	//scanner.nextLine();
	while(scanner.hasNextLine()) {
	    SAM break0 = new SAM(scanner.nextLine());
	    this.add(break0);
	}
    }

    public static void main(String[] args) {
	String fileName = "";
	if (args.length == 1) {
	    fileName = args[0];
	} else {
	    System.out.println("more/less than one argument");
	    System.exit(0);
	}
	
	SAMs breaks = new SAMs();
	breaks.run(fileName);
	//System.out.println(breaks.get(0).name);
	//System.out.println(breaks.get(0).ref);
	//String finale = breaks.printSimpleSAMs();
	//System.out.println(finale);
    }
    
}
