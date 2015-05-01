
package edu.umn.cs.kumarbio.shear;

public class Translocation extends Variant {

  private String chr2;
  private int    loc2;
  private int    size;

  // Throw null error
  public Translocation(String chr1, int loc1, String chr2, int loc2) {
    super(chr1, loc1);
    size = 0;
    this.chr2 = chr2;
    this.loc2 = loc2;
  }

  public void setSize(int size) {
    this.size = size;
  }

  @Override
  public int getSize() {
    return size;
  }

  public String getEndChr() {
    return chr2;
  }

  public int getEndLoc() {
    return loc2;
  }

}
