
package edu.umn.cs.kumarbio.shear;

public class SNP extends Variant {

  private char refBase;
  private char varBase;

  // Throw null error
  // Check for valid parameters (ACGT)
  public SNP(String chr, int loc, char refBase, char varBase) {
    super(chr, loc);
    this.refBase = refBase;
    this.varBase = varBase;
  }

  public void setRefBase(char refBase) {
    this.refBase = refBase;
  }

  public void setVarBase(char varBase) {
    this.varBase = varBase;
  }

  public char getRefBase() {
    return refBase;
  }

  public char getVarBase() {
    return varBase;
  }

  @Override
  public int getSize() {
    return 0;
  }

}
