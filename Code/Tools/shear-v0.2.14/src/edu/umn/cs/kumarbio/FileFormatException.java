/*
 * Sean Landman
 * Data Mining for Biomedical Informatics Group
 * University of Minnesota - Twin Cities
 */
package edu.umn.cs.kumarbio;

import java.io.IOException;

/**
 * Represents a file format exception.
 */
public class FileFormatException extends IOException {

  // For Serializable interface
  private static final long serialVersionUID = -2820907704538278328L;

  public FileFormatException() {
    super();
  }

  public FileFormatException(String message) {
    super(message);
  }

  public FileFormatException(Throwable cause) {
    super(cause);
  }

  public FileFormatException(String message, Throwable cause) {
    super(message, cause);
  }

}