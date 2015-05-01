
package edu.umn.cs.kumarbio.shear;

public class Inversion extends Variant {

  private int size;

  // Throw null error
  // Check for valid parameters (positive size)
  public Inversion(String chr, int loc, int size) {
    super(chr, loc);
    this.size = size;
  }

  public void setSize(int size) {
    this.size = size;
  }

  @Override
  public int getSize() {
    return size;
  }

  public int getEndLoc() {
    return (getLoc() + getSize()) - 1;
  }

}
