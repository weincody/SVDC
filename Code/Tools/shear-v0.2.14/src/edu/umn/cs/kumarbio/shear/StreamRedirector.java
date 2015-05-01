/*
 * Sean Landman
 * Data Mining for Biomedical Informatics Group
 * University of Minnesota - Twin Cities
 */

package edu.umn.cs.kumarbio.shear;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.concurrent.Callable;

/**
 * StreamRedirector is used to read the contents of an InputStream object and
 * write to an OutputStream object. This is used by the RunCommand class to run
 * system commands, and it implements Callable in order for it to be run
 * concurrently (i.e. standard error and stand output of a system command being
 * redirected at the same time). A tag can optionally be provided that will be
 * prepended to each line of the redirected content, which is useful for
 * distinguishing the system command whose output is being redirected.
 * StreamRedirector will close the InputStream object upon completion, but not
 * the OutputStream object.
 */
public class StreamRedirector implements Callable<Boolean> {

  // The InputStream object to be redirected
  public InputStream  in;

  // The OutputStream object to write to
  public OutputStream out;

  // An optional tag sting to identify the output
  public String       tag;

  /**
   * Constructs a new StreamRedirector using an InputStream object and an
   * OutputStream object. If a non-null tag string is provided then content will
   * be redirected with prepended tags.
   * 
   * @param in
   *        The InputStream object to be redirected
   * @param out
   *        The OutputStream object to write to
   * @param tag
   *        The tag string to prepend to the redirected content
   */
  public StreamRedirector(InputStream in, OutputStream out, String tag) {
    this.in = in;
    this.out = out;
    this.tag = tag;
  }

  /**
   * Constructs a new StreamRedirector using an InputStream object and an
   * OutputStream object. Content is purely redirected with no prepended tags.
   * 
   * @param in
   *        The InputStream object to be redirected
   * @param out
   *        The OutputStream object to write to
   */
  public StreamRedirector(InputStream in, OutputStream out) {
    this(in, out, null);
  }

  /**
   * Redirects the stream. If a non-null tag string is provided in this
   * invocation of StreamRedirector then content will be redirected with
   * prepended with "[SHEAR] [TAG] " on every line, where TAG is the specified
   * tag string. The InputStream object will be closed upon completion, but not
   * the OutputStream object. OutputStream can be null to ignore output.
   */
  @Override
  public Boolean call() throws IOException {

    byte[] buf = new byte[1024];
    int len;

    // For ignored output
    if (out == null) {
      while ((len = in.read(buf)) > 0) {
      }

      // For redirection with prepended tags
    } else if (tag != null) {
      boolean newLine = true;
      while ((len = in.read(buf)) > 0) {

        // Process current buffer of bytes
        for (int i = 0; i < len; i++) {

          // Write the tags if newLine = true
          if (newLine) {
            out.write(new String("[SHEAR] [" + tag + "] ").getBytes());
            newLine = false;
          }

          // Write the current byte
          out.write(buf[i]);

          // If we just wrote a "\n" set newLine = true for next time
          if (buf[i] == new String("\n").getBytes()[0]) {
            newLine = true;
          }

        }

      }

      // For regular redirection
    } else {
      while ((len = in.read(buf)) > 0) {
        out.write(buf, 0, len);
      }
    }

    // Closes the InputStream object
    in.close();

    return true;

  }
}
