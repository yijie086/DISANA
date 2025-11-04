#!/usr/bin/env groovy
// Calculate total analyzed charge with QA cuts,
// discovering runs from ANY *.hipo filenames by parsing run numbers.
// Usage:
//   groovy calcCharge.groovy runLB runUB [path]
//   groovy calcCharge.groovy [path]
// If only run limits are given, you'll be prompted for the directory.


//Load specific Modulefiles:
// 1) jdk/21.0.2   2) coatjava/13.4.0   3) qadb/3.4.0   4) groovy/4.0.20
// To run the macro, use the command:
//$COATJAVA/bin/run-groovy ./../../../../../source/rootmacros/chargeSumRuns2.groovy  path

//import org.jlab.io.hipo.HipoDataSource
import clasqa.QADB
import java.nio.file.*

String promptForPath(String prompt) {
  def console = System.console()
  if (console != null) return console.readLine(prompt)
  print(prompt); return new BufferedReader(new InputStreamReader(System.in)).readLine()
}

// Extract a run number from a filename by:
// 1) taking the LAST digit group with length >= 4 before ".hipo"
//    Examples:
//      foo_003949.hipo          -> 3949
//      DVCS_003949_02.hipo      -> 3949
//      run12_004321-file.hipo   -> 4321
//    (If multiple groups fit, we use the last one. Leading zeros ok.)
// Extract run number: take the LAST group of 4+ digits before .hipo.
// If none exists, return null (skip file).
Integer extractRunFromName(String name) {
  if (!name.toLowerCase().endsWith(".hipo")) return null
  String stem = name.substring(0, name.length() - 5) // drop ".hipo"
  def m = (stem =~ /\d{4,}/)   // only 4+ digit groups
  if (!m.find()) return null
  // collect all matches, use the last occurrence
  def groups = []
  m.reset()
  while (m.find()) groups << m.group()
  try { return groups.last().toInteger() } catch (ignored) { return null }
}


// Collect runs from a directory (non-recursive by default)
Set<Integer> collectRunsFromDir(File dir, boolean recursive = false) {
  def runs = new HashSet<Integer>()
  if (!dir?.isDirectory()) return runs
  if (!recursive) {
    dir.listFiles({ File f -> f.isFile() && f.name.endsWith(".hipo") } as FileFilter)?.each { f ->
      def rn = extractRunFromName(f.name)
      if (rn != null) runs << rn
    }
  } else {
    Files.walk(dir.toPath()).forEach { p ->
      if (Files.isRegularFile(p) && p.fileName.toString().endsWith(".hipo")) {
        def rn = extractRunFromName(p.fileName.toString())
        if (rn != null) runs << rn
      }
    }
  }
  return runs
}

void die(String msg) { System.err.println(msg); System.exit(1) }

// --------- Args ---------
Integer runLB = null, runUB = null
File dataDir = null

if (args.length == 0) {
  dataDir = new File((promptForPath("Enter directory path containing .hipo files: ") ?: "").trim())
} else if (args.length == 1) {
  dataDir = new File(args[0])
} else {
  try {
    runLB = args[0].toInteger(); runUB = args[1].toInteger()
  } catch (e) {
    die("ERROR: First two args must be integers (runLB runUB).")
  }
  if (args.length >= 3) {
    dataDir = new File(args[2])
  } else {
    dataDir = new File((promptForPath("Enter directory path containing .hipo files: ") ?: "").trim())
  }
}

if (!dataDir.exists() || !dataDir.isDirectory()) die("ERROR: '${dataDir}' is not a valid directory.")

// --------- Discover runs ---------
def runsFromDir = collectRunsFromDir(dataDir /*, recursive=false */)
if (runsFromDir.isEmpty()) die("No *.hipo files with parseable run numbers found in ${dataDir}.")

if (runLB == null || runUB == null) {
  runLB = runsFromDir.min(); runUB = runsFromDir.max()
  println "Inferred run range from files: ${runLB} to ${runUB}"
} else if (runLB > runUB) {
  die("ERROR: runLB (${runLB}) must be <= runUB (${runUB}).")
}

println "Scanning: ${dataDir.absolutePath}  |  Found ${runsFromDir.size()} run IDs in filenames."
def allowedRuns = runsFromDir.findAll { it >= runLB && it <= runUB } as Set<Integer>
if (allowedRuns.isEmpty()) die("No runs in directory fall within ${runLB}..${runUB}.")

// --------- QA setup ---------
println "run range: $runLB to $runUB (restricted to runs present in directory)"
QADB qa = new QADB("latest", runLB, runUB)

// Enable desired defect checks
qa.checkForDefect("TotalOutlier");
qa.checkForDefect("TerminalOutlier");
qa.checkForDefect("MarginalOutlier");
qa.checkForDefect("SectorLoss");
// qa.checkForDefect("LowLiveTime");
qa.checkForDefect("Misc");
qa.checkForDefect("TotalOutlierFT");
qa.checkForDefect("TerminalOutlierFT");
qa.checkForDefect("MarginalOutlierFT");
qa.checkForDefect("LossFT");
qa.checkForDefect("BSAWrong");
qa.checkForDefect("BSAUnknown");
qa.checkForDefect("TSAWrong");
qa.checkForDefect("TSAUnknown");
qa.checkForDefect("DSAWrong");
qa.checkForDefect("DSAUnknown");
qa.checkForDefect("ChargeHigh");
qa.checkForDefect("ChargeNegative");
qa.checkForDefect("ChargeUnknown");
qa.checkForDefect("PossiblyNoBeam");

// --------- Iterate only over runs that exist in dir ---------
int runnum, evnum
qa.getQaTree().each { runnumStr, runTree ->
  runnum = runnumStr.toInteger()
  if (!allowedRuns.contains(runnum)) return

  double beforeTotal = qa.getAccumulatedCharge()
  double beforePos   = qa.getAccumulatedChargeHL(1)
  double beforeNeg   = qa.getAccumulatedChargeHL(-1)
  double before0     = qa.getAccumulatedChargeHL(0)

  runTree.each { filenumStr, fileTree ->
    evnum = fileTree['evnumMin']
    if (qa.pass(runnum, evnum)) {
      qa.accumulateCharge()
      qa.accumulateChargeHL()
    }
  }

  double runTotal = qa.getAccumulatedCharge() - beforeTotal
  double runPos   = qa.getAccumulatedChargeHL(1) - beforePos
  double runNeg   = qa.getAccumulatedChargeHL(-1) - beforeNeg
  double run0     = qa.getAccumulatedChargeHL(0) - before0

  println "${runnum},${runTotal},${runPos},${runNeg},${run0}"
}

// --------- Totals ---------
println "\ntotal accumulated charge: " + qa.getAccumulatedCharge() + " (nC)"
println "total no hel accumulated charge: " + qa.getAccumulatedChargeHL(0) + " (nC)"
println "total pos hel accumulated charge: " + qa.getAccumulatedChargeHL(1) + " (nC)"
println "total neg hel accumulated charge: " + qa.getAccumulatedChargeHL(-1) + " (nC)"
println "total unassigned hel accumulated charge: " + qa.getAccumulatedChargeHL(0) + " (nC)"
