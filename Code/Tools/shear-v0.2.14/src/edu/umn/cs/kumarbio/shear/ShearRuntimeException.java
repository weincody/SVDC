/*
 * Sean Landman
 * Data Mining for Biomedical Informatics Group
 * University of Minnesota - Twin Cities
 */

package edu.umn.cs.kumarbio.shear;

/**
 * Represents a runtime exception (i.e. unexpected behavior) for SHEAR
 */
public class ShearRuntimeException extends RuntimeException {

  // For Serializable interface
  private static final long serialVersionUID = 3042370574173739731L;

  public ShearRuntimeException() {
    super();
  }

  public ShearRuntimeException(String message) {
    super(message);
  }

  public ShearRuntimeException(Throwable cause) {
    super(cause);
  }

  public ShearRuntimeException(String message, Throwable cause) {
    super(message, cause);
  }

}
