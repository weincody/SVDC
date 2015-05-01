/*
 * Sean Landman Data Mining for Biomedical Informatics Group University of
 * Minnesota - Twin Cities
 */

package edu.umn.cs.kumarbio.shear;

import java.io.File;
import java.io.IOException;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.FutureTask;

/**
 * RunPicardCommand is used to run a PicardTools command. It contains only the
 * static run method, which takes a string specifying the PicardTools JAR file
 * to use and runs it, redirecting the output and error streams to be prepended
 * with a tag and sent to standard output.
 */
public class RunPicardCommand {

  /**
   * Runs the designated PicardTools program.
   * 
   * @param picardDir
   *          String specifying the directory containing the PicardTools JAR
   *          files
   * @param picardProgram
   *          String specifying the PicardTools JAR file to run
   * @param picardArgs
   *          Array of arguments for the PicardTools command
   * @param debug
   *          Whether or not debugging mode is on
   */
  public static void run(String picardDir, String picardProgram, String[] picardArgs, boolean debug) {

    // Print out the command being executed
    if (debug) {
      System.out.println("[SHEAR] ------------------------------------------------------------------------");
      System.out.print("[SHEAR] PICARD-COMMAND: " + picardProgram);
      for (int i = 0; i < picardArgs.length; i++) {
        System.out.print(" " + picardArgs[i]);
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
      FutureTask<Boolean> redirector =
          new FutureTask<Boolean>(new StreamRedirector(pipedIn, redirectionTarget, "PICARD"));
      new Thread(redirector).start();

      // Verify correct picardProgram argument
      if (!(picardProgram.equals("AddOrReplaceReadGroups") || picardProgram.equals("CleanSam") || picardProgram
          .equals("FixMateInformation"))) {
        throw new ShearRuntimeException("Illegal Picard tool passed to RunPicardCommand.run method: " + picardProgram);
      }

      // Get JAR path
      String picardJarPath = picardDir;
      if (new File(picardDir + "/picard.jar").exists()) {
        picardJarPath += "/" + "picard.jar";
      } else {
        picardJarPath += "/" + picardProgram + ".jar";
      }

      // Run PicardTools command
      URLClassLoader picardClassLoader =
          new URLClassLoader(new URL[] { new URL("file", null, new File(picardJarPath).getAbsolutePath()) }, null);
      Class<?> picardClass;
      try {
        picardClass = picardClassLoader.loadClass("picard.sam." + picardProgram);
      } catch (ClassNotFoundException e) {
        picardClass = picardClassLoader.loadClass("net.sf.picard.sam." + picardProgram);
      }
      Method mainMethod = picardClass.getMethod("instanceMain", new Class<?>[] { String[].class });
      mainMethod.invoke(picardClass.newInstance(), new Object[] { picardArgs });
      picardClassLoader.close();

      // Closed tagged out/err streams
      System.out.close();
      System.err.close();

      // Wait for the streams
      redirector.get();

      // Print ending separator
      if (debug) {
        System.out.println("[SHEAR] ------------------------------------------------------------------------");
      }

      // Catch IO and invocation errors
    } catch (InterruptedException e) {
      throw new ShearRuntimeException("Current thread was interrupted while waiting for stream redirection", e);
    } catch (ExecutionException e) {
      throw new ShearRuntimeException("Stream redirection thread encountered an error", e);
    } catch (IOException e) {
      throw new ShearRuntimeException(
          "Error encountered with PipedOutputStream and/or PipedInputStream objects for stream redirection",
          e);
    } catch (ClassNotFoundException e) {
      throw new ShearRuntimeException("Unable to load necessary Picard class", e);
    } catch (NoSuchMethodException e) {
      throw new ShearRuntimeException("Unable to load necessary Picard method", e);
    } catch (IllegalAccessException e) {
      throw new ShearRuntimeException("Error encountered with accesing Picard method", e);
    } catch (InstantiationException e) {
      throw new ShearRuntimeException("Error encountered with instantiating Picard object", e);
    } catch (InvocationTargetException e) {
      throw new ShearRuntimeException("Error encountered with invoking Picard method", e);

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
