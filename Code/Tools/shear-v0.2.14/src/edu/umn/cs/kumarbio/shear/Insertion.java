
package edu.umn.cs.kumarbio.shear;

public class Insertion extends Variant {

  private String insertedBases;

  // Throw null error
  // Check for valid parameters (positive size)
  public Insertion(String chr, int loc, String insertedBases) {
    super(chr, loc);
    this.insertedBases = insertedBases;
  }

  public void setInsertedBases(String insertedBases) {
    this.insertedBases = insertedBases;
  }

  public String getInsertedBases() {
    return insertedBases;
  }

  @Override
  public int getSize() {
    return getInsertedBases().length();
  }

}
