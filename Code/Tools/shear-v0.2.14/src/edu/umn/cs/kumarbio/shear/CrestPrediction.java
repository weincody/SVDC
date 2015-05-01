
package edu.umn.cs.kumarbio.shear;

public class CrestPrediction {

  public String type;
  public String bp1region;
  public int    bp1;
  public String bp1dir;
  public int    bp1sc;
  public int    bp1span;
  public String bp2region;
  public int    bp2;
  public String bp2dir;
  public int    bp2sc;
  public int    bp2span;
  public double heterogeneityPercent;

  public CrestPrediction(String type, String bp1region, int bp1, String bp1dir, int bp1sc, int bp1span,
      String bp2region, int bp2, String bp2dir, int bp2sc, int bp2span, double heterogeneityPercent) {
    this.type = type;
    this.bp1region = bp1region;
    this.bp1 = bp1;
    this.bp1dir = bp1dir;
    this.bp1sc = bp1sc;
    this.bp1span = bp1span;
    this.bp2region = bp2region;
    this.bp2 = bp2;
    this.bp2dir = bp2dir;
    this.bp2sc = bp2sc;
    this.bp2span = bp2span;
    this.heterogeneityPercent = heterogeneityPercent;
  }
}
