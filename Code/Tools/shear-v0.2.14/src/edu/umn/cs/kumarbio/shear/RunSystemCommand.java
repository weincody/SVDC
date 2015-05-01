/*
 * Sean Landman
 * Data Mining for Biomedical Informatics Group
 * University of Minnesota - Twin Cities
 */

package edu.umn.cs.kumarbio.shear;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.util.List;
import java.util.concurrent.FutureTask;

/**
 * RunSystemCommand is used to execute a system command and direct the output.
 * It contains only the static run method, which takes a ProcessBuilder object
 * and runs it, possibly redirecting the output and error streams to a specified
 * OutputStream object, or returning output in a string. Any stream can also
 * have each line prepended with a tag to indicate the external program being
 * run.
 */
public class RunSystemCommand {

  /**
   * Runs the specified process and returns a string containing the output and
   * error stream results.
   * 
   * @param pb
   *        The ProcessBuilder object containing the command to run
   * @param debug
   *        Whether or not debugging mode is on
   * @return The output and error stream results of the process
   */
  public static String run(ProcessBuilder pb, boolean debug) {
    try {
      pb.redirectErrorStream(true);
      ByteArrayOutputStream result = new ByteArrayOutputStream();
      Process p = pb.start();
      FutureTask<Boolean> redirector = new FutureTask<Boolean>(new StreamRedirector(p.getInputStream(), result));
      new Thread(redirector).start();
      redirector.get();
      p.waitFor();
      return result.toString();
    } catch (Exception e) {
      throw new ShearRuntimeException(e);
    }
  }

/*  *//**
   * Runs the specified process and redirects both the output and error
   * streams
   * to the specified OutputStream object.
   * 
   * @param pb
   *        The ProcessBuilder object containing the command to run
   * @param combined
   *        The OutputStream to redirect both the output and error streams to
   * @param debug
   *        Whether or not debugging mode is on
   */
  /*
   * public static void run(ProcessBuilder pb, OutputStream combined, boolean
   * debug) {
   * run(pb, combined, null, null, null, debug);
   * }
   *//**
   * Runs the specified process and redirects both the output and error
   * streams
   * to the specified OutputStream object.
   * 
   * @param pb
   *        The ProcessBuilder object containing the command to run
   * @param combined
   *        The OutputStream to redirect both the output and error streams to
   * @param combinedTag
   *        The tag string to prepend to the redirected content
   * @param debug
   *        Whether or not debugging mode is on
   */
  /*
   * public static void run(ProcessBuilder pb, OutputStream combined, String
   * combinedTag, boolean debug) {
   * run(pb, combined, combinedTag, null, null, debug);
   * }
   *//**
   * Runs the specified process and redirects both the output and error
   * streams
   * to their respective specified OutputStream objects.
   * 
   * @param pb
   *        The ProcessBuilder object containing the command to run
   * @param stout
   *        The OutputStream to redirect the output stream to
   * @param sterr
   *        The OutputStream to redirect the error stream to
   * @param debug
   *        Whether or not debugging mode is on
   */
  /*
   * public static void run(ProcessBuilder pb, OutputStream stdout, OutputStream
   * stderr, boolean debug) {
   * run(pb, stdout, null, stderr, null, debug);
   * }
   *//**
   * Runs the specified process and redirects both the output and error
   * streams
   * to their respective specified OutputStream objects.
   * 
   * @param pb
   *        The ProcessBuilder object containing the command to run
   * @param stout
   *        The OutputStream to redirect the output stream to
   * @param outTag
   *        The tag string to prepend to the redirected output stream content
   * @param sterr
   *        The OutputStream to redirect the error stream to
   * @param debug
   *        Whether or not debugging mode is on
   */
  /*
   * public static void run(ProcessBuilder pb, OutputStream stdout, String
   * outTag, OutputStream stderr, boolean debug) {
   * run(pb, stdout, outTag, stderr, null, debug);
   * }
   *//**
   * Runs the specified process and redirects both the output and error
   * streams
   * to their respective specified OutputStream objects.
   * 
   * @param pb
   *        The ProcessBuilder object containing the command to run
   * @param stout
   *        The OutputStream to redirect the output stream to
   * @param sterr
   *        The OutputStream to redirect the error stream to
   * @param errTag
   *        The tag string to prepend to the redirected error stream content
   * @param debug
   *        Whether or not debugging mode is on
   */
  /*
   * public static void run(ProcessBuilder pb, OutputStream stdout, OutputStream
   * stderr, String errTag, boolean debug) {
   * run(pb, stdout, null, stderr, errTag, debug);
   * }
   */

  /**
   * Runs the specified process and redirects both the output and error streams
   * to their respective specified OutputStream objects.
   * 
   * @param pb
   *        The ProcessBuilder object containing the command to run
   * @param stout
   *        The OutputStream to redirect the output stream to
   * @param outTag
   *        The tag string to prepend to the redirected output stream content
   * @param sterr
   *        The OutputStream to redirect the error stream to
   * @param errTag
   *        The tag string to prepend to the redirected error stream content
   * @param debug
   *        Whether or not debugging mode is on
   */
  public static void run(ProcessBuilder pb, OutputStream stdout, String outTag, OutputStream stderr, String errTag,
      boolean debug) {

    try {

      // Redirect error stream to output stream if the same OutputStream object is given for them
      if (stderr == stdout) {
        pb.redirectErrorStream(true);
      }

      // Print out the command being executed
      if (debug) {
        System.out.println("[SHEAR] ------------------------------------------------------------------------");
        List<String> tokens = pb.command();
        System.out.print("[SHEAR] SYSTEM-COMMAND:");
        for (int i = 0; i < tokens.size(); i++) {
          System.out.print(" " + tokens.get(i));
        }
        System.out.println();
        System.out.println("[SHEAR] ------------------------------------------------------------------------");
      }

      // Start the process
      Process p = pb.start();

      // Redirect the output stream
      FutureTask<Boolean> redirector = new FutureTask<Boolean>(new StreamRedirector(p.getInputStream(), stdout, outTag));
      new Thread(redirector).start();

      // Redirect the error stream, and wait for it
      if (!pb.redirectErrorStream()) {
        FutureTask<Boolean> errRedirector = new FutureTask<Boolean>(new StreamRedirector(p.getErrorStream(), stderr,
            errTag));
        new Thread(errRedirector).start();
        errRedirector.get();
      }

      // Wait for the output stream and for the process to finish
      redirector.get();
      p.waitFor();

      // Print ending separator
      if (debug) {
        System.out.println("[SHEAR] ------------------------------------------------------------------------");
      }

    } catch (Exception e) {
      throw new ShearRuntimeException(e);
    }

  }

}
