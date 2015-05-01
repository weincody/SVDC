package edu.umn.cs.kumarbio.shear;

public class SdiLine {
  public String region;
  public int    loc;
  public String origBases;
  public String newBases;

  public SdiLine(String region, int loc, String origBases, String newBases) {
    this.region = region;
    this.loc = loc;
    this.origBases = origBases;
    this.newBases = newBases;
  }
}
