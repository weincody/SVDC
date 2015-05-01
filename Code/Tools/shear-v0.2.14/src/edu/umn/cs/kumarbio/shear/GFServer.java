/*
 * Sean Landman Data Mining for Biomedical Informatics Group University of
 * Minnesota - Twin Cities
 */

package edu.umn.cs.kumarbio.shear;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Scanner;
import java.util.concurrent.ExecutionException;

public class GFServer extends Thread {

  private final String         prefix;
  private final String         twoBitFile;
  private int                  portNum;
  private final ProcessBuilder pb;
  private Process              proc;
  private boolean              debug;
  private boolean              validTwoBitFile;

  // TO DO: Need to fail nice for non-unix by checking for netstat/pgrep
  public GFServer(String prefix, String twoBitFile, boolean debug)
      throws InterruptedException, IOException, ExecutionException {
    this(prefix, twoBitFile, debug, 49152);
  }

  // TO DO: Need to fail nice for non-unix by checking for netstat/pgrep
  public GFServer(String prefix, String twoBitFile, boolean debug, int startingPortNum)
      throws InterruptedException, IOException, ExecutionException {
    this.prefix = prefix;
    this.twoBitFile = twoBitFile;
    this.debug = debug;
    validTwoBitFile = true;
    clearLogFiles();
    portNum = getNextAvailablePort(startingPortNum);
    if (portNum == -1) {
      throw new InterruptedException("Could not locate a port number to use!");
    }
    pb =
        new ProcessBuilder("gfServer", "start", "localhost", "" + portNum, twoBitFile, "-canStop", "-log=" + prefix
            + ".gfServer.log");
  }

  private int getNextAvailablePort(int startingPortNum) throws InterruptedException, IOException, ExecutionException {
    // String portCheck = RunCommand.run(new ProcessBuilder("lsof", "-i",
    // "@0.0.0.0:" + port));
    // while (portCheck.length() != 0) {
    // port++;
    // portCheck = RunCommand.run(new ProcessBuilder("lsof", "-i", "@0.0.0.0:" +
    // port));
    // }
    String processCheck = RunSystemCommand.run(new ProcessBuilder("pgrep", "-l", "-f", "gfServer"), debug);
    String portCheck = RunSystemCommand.run(new ProcessBuilder("netstat", "-an"), debug);
    for (int i = startingPortNum; i <= 65535; i++) {
      // String portCheck = RunSystemCommand.run(new ProcessBuilder("lsof",
      // "-i", "@0.0.0.0:" + i), debug);
      if (portCheck.contains("" + i) || processCheck.contains("" + i)) {
        continue;
      } else {
        return i;
      }
    }
    return -1;
    // for (int i = 49152; i <= 65535; i++) {
    // ServerSocket serverSocket = null;
    // DatagramSocket datagramSocket = null;
    // try {
    // serverSocket = new ServerSocket(i);
    // serverSocket.setReuseAddress(true);
    // datagramSocket = new DatagramSocket(i);
    // datagramSocket.setReuseAddress(true);
    // return i;
    // } catch (IOException e) {
    // continue;
    // } finally {
    // if (datagramSocket != null) {
    // datagramSocket.close();
    // }
    // if (serverSocket != null) {
    // try {
    // serverSocket.close();
    // } catch (IOException e) {
    // // Can never happen
    // throw new InterruptedException("Could not locate a port number to use!");
    // }
    // }
    // }
    // }
    // return -1;
  }

  // private boolean isPortFree(int port) throws InterruptedException,
  // IOException, ExecutionException {
  // String portCheck = RunCommand.run(new ProcessBuilder("pgrep", "-l", "-f",
  // "gfServer"));
  // if (portCheck.contains("" + port)) {
  // return false;
  // } else {
  // return true;
  // }
  // }

  @Override
  public void run() {
    try {
      pb.redirectErrorStream(true);
      proc = pb.start();
      InputStream in = proc.getInputStream();
      FileOutputStream out = new FileOutputStream(prefix + ".gfServer.out");
      byte[] buf = new byte[1024];
      int len;
      while ((len = in.read(buf)) > 0) {
        out.write(buf, 0, len);
        out.flush();
      }
      in.close();
      out.close();
      proc.waitFor();
    } catch (InterruptedException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  // Returns false if a socket/connection error is encountered in the log file
  public boolean waitForReady() throws InterruptedException, IOException, ExecutionException {
    while (proc == null) {
      sleep(1000);
    }
    boolean logFileReady = false;
    String serverCheck =
        RunSystemCommand.run(new ProcessBuilder("gfServer", "files", "localhost", "" + portNum), debug);
    while (!logFileReady || !serverCheck.trim().equals(twoBitFile)) {
      sleep(1000);
      Scanner sc = new Scanner(new File(prefix + ".gfServer.log"));
      int ctr = 0;
      while (ctr < 10 && sc.hasNextLine()) {
        ctr++;
        String line = sc.nextLine();
        if (line.contains("Couldn't bind socket to") || line.contains("Error accepting the connection")) {
          sc.close();
          return false;
        } else if (line.contains("doesn't have a valid twoBitSig") || line.contains("Unrecognized file type")) {
          validTwoBitFile = false;
          sc.close();
          return false;
        } else if (line.contains("Server ready for queries")) {
          logFileReady = true;
        }
      }
      sc.close();
      serverCheck = RunSystemCommand.run(new ProcessBuilder("gfServer", "files", "localhost", "" + portNum), debug);
    }
    return true;
  }

  public int getPortNum() {
    return portNum;
  }

  public String getTwoBitFile() {
    return twoBitFile;
  }

  public boolean isValidTwoBitFile() {
    return validTwoBitFile;
  }

  public void close() throws InterruptedException {
    if (proc != null) {
      proc.destroy();
    }
  }

  public void clearLogFiles() {
    new File(prefix + ".gfServer.log").delete();
    new File(prefix + ".gfServer.out").delete();
  }

}
