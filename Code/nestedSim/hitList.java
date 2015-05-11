/*
  Code for comparing corresponding Socrates outputs
*/

import java.util.*;
import java.io.*;
import java.math.*;
import java.lang.*;


public class hitList {
    public ArrayList<hitTest> book;

    public hitList(){
	book = new ArrayList<hitTest>();
    }

    public void add(hitTest break0){
	boolean added = false;
	if(book.size() == 0){
	    book.add(break0);
	    added = true;
	    return;
	}
	book.add(break0);
    }

    public hitTest get(int i){
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
	String temp = "";
	while(scanner.hasNextLine()) {
	    temp = scanner.nextLine();
	    hitTest hitTest0 = new hitTest(temp);
	    this.add(hitTest0);
	}
	//this.printBook();		
    }

    public static void main(String[] args) {
	String fileName = "";
	if (args.length == 1) {
	    fileName = args[0];
	} else {
	    System.out.println("more than one argument");
	    System.exit(0);
	}
	
	hitList hitList0 = new hitList();
	hitList0.run(fileName);
    }
    
}
