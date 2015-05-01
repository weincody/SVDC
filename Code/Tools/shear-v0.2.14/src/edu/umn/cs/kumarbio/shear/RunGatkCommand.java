/*
 * Sean Landman
 * Data Mining for Biomedical Informatics Group
 * University of Minnesota - Twin Cities
 */

package edu.umn.cs.kumarbio.shear;

import java.io.IOException;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.FutureTask;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;

/**
 * RunGatkCommand is used to run a GATK command. It contains only the static run
 * method, which takes a string array specifying the GATK parameters and runs
 * it, redirecting the output and error streams to be prepended with a tag and
 * sent to standard output.
 */
public class RunGatkCommand {

  /**
   * Runs the designated GATK command.
   * 
   * @param gatkArgs
   *        Array of arguments for the GATK command
   * @param debug
   *        Whether or not debugging mode is on
   */
  public static void run(String[] gatkArgs, boolean debug) {

    // Print out the command being executed
    if (debug) {
      System.out.println("[SHEAR] ------------------------------------------------------------------------");
      System.out.print("[SHEAR] GATK-COMMAND:");
      for (int i = 0; i < gatkArgs.length; i++) {
        System.out.print(" " + gatkArgs[i]);
      }
      System.out.println();
      System.out.println("[SHEAR] ------------------------------------------------------------------------");
    }

    // Save original standard out/err streams
    PrintStream origSystemOut = System.out;
    PrintStream origSystemErr = System.err;

    // Set up piped output->input streams
    PipedOutputStream pipedOut = new PipedOutputStream();
    PipedInputStream pipedIn = new PipedInputStream();

    // Wrap piped output in a PrintStream object
    PrintStream printStream = new PrintStream(pipedOut);

    try {

      // Connect piped output->input streams
      pipedOut.connect(pipedIn);

      // Point tagged out/err to piped output
      System.setOut(printStream);
      System.setErr(printStream);

      // Redirect the output stream in separate thread
      PrintStream redirectionTarget = debug ? origSystemOut : null;
      FutureTask<Boolean> redirector = new FutureTask<Boolean>(new StreamRedirector(pipedIn, redirectionTarget, "GATK"));
      new Thread(redirector).start();

      // Run GATK command
      CommandLineProgram.start(new CommandLineGATK(), gatkArgs);

      // Closed tagged out/err streams
      System.out.close();
      System.err.close();

      // Wait for the streams
      redirector.get();

      // Print ending separator
      if (debug) {
        System.out.println("[SHEAR] ------------------------------------------------------------------------");
      }

      // Catch IO errors
    } catch (InterruptedException e) {
      throw new ShearRuntimeException("Current thread was interrupted while waiting for stream redirection", e);
    } catch (ExecutionException e) {
      throw new ShearRuntimeException("Stream redirection thread encountered an error", e);
    } catch (IOException e) {
      throw new ShearRuntimeException(
          "Error encountered with PipedOutputStream and/or PipedInputStream objects for stream redirection", e);
    } catch (Exception e) {
      throw new ShearRuntimeException("Error with running GATK command", e);

    } finally {

      // Close tagged out/err streams
      System.out.close();
      System.err.close();

      // Closed piped output->input streams
      printStream.close();

      // Reset original standard out/err streams
      System.setOut(origSystemOut);
      System.setErr(origSystemErr);

    }

  }

}
