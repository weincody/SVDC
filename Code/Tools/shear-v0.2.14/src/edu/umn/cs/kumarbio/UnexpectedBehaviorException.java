/*
 * Sean Landman
 * Data Mining for Biomedical Informatics Group
 * University of Minnesota - Twin Cities
 */

package edu.umn.cs.kumarbio;

/**
 * Represents an unexpected behavior exception.
 */
public class UnexpectedBehaviorException extends RuntimeException {

  // For Serializable interface
  private static final long serialVersionUID = -1148608975110084758L;

  public UnexpectedBehaviorException() {
    super();
  }

  public UnexpectedBehaviorException(String message) {
    super(message);
  }

  public UnexpectedBehaviorException(Throwable cause) {
    super(cause);
  }

  public UnexpectedBehaviorException(String message, Throwable cause) {
    super(message, cause);
  }

}
