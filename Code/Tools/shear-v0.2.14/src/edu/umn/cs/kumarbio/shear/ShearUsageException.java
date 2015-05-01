/*
 * Sean Landman
 * Data Mining for Biomedical Informatics Group
 * University of Minnesota - Twin Cities
 */

package edu.umn.cs.kumarbio.shear;

/**
 * Represents a runtime exception (i.e. unexpected behavior) for SHEAR
 */
public class ShearUsageException extends Exception {

  // For Serializable interface
  private static final long serialVersionUID = 2107953043229978796L;

  public ShearUsageException() {
    super();
  }

  public ShearUsageException(String message) {
    super(message);
  }

  public ShearUsageException(Throwable cause) {
    super(cause);
  }

  public ShearUsageException(String message, Throwable cause) {
    super(message, cause);
  }

}
