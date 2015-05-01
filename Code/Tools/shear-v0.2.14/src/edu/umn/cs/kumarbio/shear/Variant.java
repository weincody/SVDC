
package edu.umn.cs.kumarbio.shear;

public abstract class Variant implements Comparable<Variant> {

  private String  chr;
  private int     loc;
  private double  heterogeneity;
  private boolean isIndel;

  // Throw null error
  public Variant(String chr, int loc) {
    this.chr = chr;
    this.loc = loc;
    isIndel = false;
  }

  public void setChr(String chr) {
    this.chr = chr;
  }

  public void setLoc(int loc) {
    this.loc = loc;
  }

  public void setIndel(boolean isIndel) {
    this.isIndel = isIndel;
  }

  public void setHeterogeneity(double heterogeneity) {
    this.heterogeneity = heterogeneity;
  }

  public String getChr() {
    return chr;
  }

  public int getLoc() {
    return loc;
  }

  // Throw error for null heterogeneity?
  public double getHeterogeneity() {
    return heterogeneity;
  }

  public boolean isIndel() {
    return isIndel;
  }

  public abstract int getSize();

  public boolean overlapsWith(Variant v) {
    int v1start = getLoc();
    int v1end = (getLoc() + getSize()) - 1;
    int v2start = v.getLoc();
    int v2end = (v.getLoc() + v.getSize()) - 1;
    return ((v1start >= v2start) && (v1start <= v2end)) || ((v1end >= v2start) && (v1end <= v2end));
  }

  public boolean adjacentTo(Variant v) {
    int v1start = getLoc();
    int v1end = (getLoc() + getSize()) - 1;
    int v2start = v.getLoc();
    int v2end = (v.getLoc() + v.getSize()) - 1;
    return ((v1start >= (v2start - 1)) && (v1start <= (v2end + 1)))
        || ((v1end >= (v2start - 1)) && (v1end <= (v2end + 1)));
  }

  @Override
  public int compareTo(Variant v) {
    // Need to improve! Check for equality vs. just same loc
    int comparison = getLoc() - v.getLoc();
    return comparison == 0 ? -1 : comparison;
  }

}
