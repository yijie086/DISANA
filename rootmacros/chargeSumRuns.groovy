#!/usr/bin/env groovy
// Calculate total analyzed charge with QA cuts.
// Modes:
//   1) run-only:      runLB runUB
//   2) path-only:     [path]
//   3) combined:      runLB runUB [path]
// Optional flags anywhere: --recursive  --exclude "3211,3215-3220"  --exclude-file file
//
// Typical modules on ifarm: jdk/21.0.2, coatjava/13.4.0, qadb/3.4.0, groovy/4.0.20
// Run examples:
//   $COATJAVA/bin/run-groovy chargeSumRuns.groovy 5681 5757
//   $COATJAVA/bin/run-groovy chargeSumRuns.groovy /path/to/DVCSWagon/
//   $COATJAVA/bin/run-groovy chargeSumRuns.groovy 5681 5757 /path/to/DVCSWagon/ --exclude "5690,5705-5708"

import clasqa.QADB
import java.nio.file.*

// ---------- helpers ----------
String promptForPath(String prompt) {
  def c = System.console()
  if (c != null) return c.readLine(prompt)
  print(prompt); return new BufferedReader(new InputStreamReader(System.in)).readLine()
}

boolean isInt(String s) {
  try { s?.toInteger(); return true } catch (ignored) { return false }
}

// Extract run number: LAST group of 4+ digits before .hipo; else null.
Integer extractRunFromName(String name) {
  if (!name?.toLowerCase()?.endsWith(".hipo")) return null
  String stem = name.substring(0, name.length() - 5)
  def m = (stem =~ /\d{4,}/)
  if (!m.find()) return null
  def gs = []
  m.reset(); while (m.find()) gs << m.group()
  try { return gs.last().toInteger() } catch (ignored) { return null }
}

// Collect runs from a directory (case-insensitive *.hipo). Set recursive=true to walk subdirs.
Set<Integer> collectRunsFromDir(File dir, boolean recursive = false) {
  def runs = new LinkedHashSet<Integer>()
  if (!dir?.isDirectory()) return runs
  if (!recursive) {
    dir.eachFileMatch(~/(?i).*\.hipo$/) { f ->
      if (f.isFile()) {
        def rn = extractRunFromName(f.name)
        if (rn != null) runs << rn
      }
    }
  } else {
    dir.eachFileRecurse { f ->
      if (f.isFile() && f.name =~ /(?i).*\.hipo$/) {
        def rn = extractRunFromName(f.name)
        if (rn != null) runs << rn
      }
    }
  }
  return runs
}

void die(String msg) { System.err.println(msg); System.exit(1) }

// Parse "3211,3215-3220,3987"
Set<Integer> parseRunList(String spec) {
  def out = new LinkedHashSet<Integer>()
  if (!spec) return out
  spec.split(/[,\s]+/).each { tok ->
    if (!tok) return
    if (tok.contains("-")) {
      def parts = tok.split("-", 2)
      try {
        int a = parts[0].trim().toInteger()
        int b = parts[1].trim().toInteger()
        int lo = Math.min(a,b), hi = Math.max(a,b)
        (lo..hi).each { out << it }
      } catch (ignored) {}
    } else {
      try { out << tok.trim().toInteger() } catch (ignored) {}
    }
  }
  return out
}

Set<Integer> loadExcludeFile(String path) {
  def out = new LinkedHashSet<Integer>()
  if (!path) return out
  def f = new File(path)
  if (!f.exists() || !f.isFile()) die("ERROR: exclude file '${path}' not found.")
  f.eachLine { line -> out.addAll(parseRunList(line)) }
  return out
}

// ---------- parse flags (can appear anywhere) ----------
def argList = (args as List).collect { it as String }
boolean recursive = false
String excludeSpec = null
String excludeFile = null

for (int i = 0; i < argList.size(); ) {
  def a = argList[i]
  if (a == "--recursive") { recursive = true; argList.remove(i) }
  else if (a == "--exclude" && i+1 < argList.size()) { excludeSpec = argList[i+1]; argList.remove(i+1); argList.remove(i) }
  else if (a == "--exclude-file" && i+1 < argList.size()) { excludeFile = argList[i+1]; argList.remove(i+1); argList.remove(i) }
  else { i++ }
}

// ---------- decide mode from remaining positionals ----------
Integer runLB = null, runUB = null
File dataDir = null
boolean runOnlyMode = false
boolean pathOnlyMode = false
boolean combinedMode = false

if (argList.size() == 0) {
  // no args -> ask for path (path-only)
  dataDir = new File((promptForPath("Enter directory path containing .hipo files: ") ?: "").trim())
  pathOnlyMode = true
} else if (argList.size() == 1) {
  // single positional: either a path OR a run number (we treat as path)
  if (isInt(argList[0])) {
    die("ERROR: Received one integer ('${argList[0]}'). For run-only mode provide: runLB runUB")
  }
  dataDir = new File(argList[0])
  pathOnlyMode = true
} else if (argList.size() >= 2 && isInt(argList[0]) && isInt(argList[1])) {
  // have runLB runUB
  runLB = argList[0].toInteger()
  runUB = argList[1].toInteger()
  if (runLB > runUB) die("ERROR: runLB (${runLB}) must be <= runUB (${runUB}).")
  if (argList.size() >= 3) {
    dataDir = new File(argList[2])
    combinedMode = true
  } else {
    runOnlyMode = true
  }
} else {
  die("ERROR: Unrecognized arguments. Use either [runLB runUB] or [path] or [runLB runUB path].")
}

// ---------- build allowedRuns ----------
Set<Integer> allowedRuns
if (runOnlyMode) {
  // No filesystem: use the entire range
  allowedRuns = (runLB..runUB).toSet()
  println "Mode: run-only. Using runs ${runLB}..${runUB} (no directory filtering)."
} else {
  if (!dataDir.exists() || !dataDir.isDirectory()) die("ERROR: '${dataDir}' is not a valid directory.")
  def runsFromDir = collectRunsFromDir(dataDir, recursive)
  if (runsFromDir.isEmpty()) die("No *.hipo files with parseable run numbers found in ${dataDir}.")

  if (pathOnlyMode) {
    runLB = runsFromDir.min(); runUB = runsFromDir.max()
    println "Mode: path-only. Inferred run range from files: ${runLB} to ${runUB}"
    allowedRuns = runsFromDir
  } else {
    println "Mode: combined. Directory has ${runsFromDir.size()} run IDs."
    allowedRuns = runsFromDir.findAll { it >= runLB && it <= runUB } as Set<Integer>
    if (allowedRuns.isEmpty()) die("No runs in directory fall within ${runLB}..${runUB}.")
  }
  println "Scanning: ${dataDir.absolutePath}"
}

// ---------- apply exclusions (all modes) ----------
def excludeRuns = new LinkedHashSet<Integer>()
excludeRuns.addAll(parseRunList(excludeSpec))
excludeRuns.addAll(loadExcludeFile(excludeFile))
excludeRuns.addAll(parseRunList(System.getenv("EXCLUDE_RUNS")))
if (!excludeRuns.isEmpty()) {
  def before = allowedRuns.size()
  allowedRuns.removeAll(excludeRuns)
  println "Excluded ${excludeRuns.size()} specified runs; kept ${allowedRuns.size()} of ${before}."
  if (allowedRuns.isEmpty()) die("No runs left to analyze after exclusions.")
}

// ---------- QA setup ----------
println "run range considered: $runLB to $runUB"
println "number of runs to analyze: ${allowedRuns.size()}"
QADB qa = new QADB("latest", runLB, runUB)

// Enable desired defect checks
qa.checkForDefect("TotalOutlier");
qa.checkForDefect("TerminalOutlier");
qa.checkForDefect("MarginalOutlier");
qa.checkForDefect("SectorLoss");
//qa.checkForDefect("LowLiveTime"); // enabled per your note for RGA spring 2018
qa.checkForDefect("Misc");
qa.checkForDefect("TotalOutlierFT");
qa.checkForDefect("TerminalOutlierFT");
qa.checkForDefect("MarginalOutlierFT");
qa.checkForDefect("LossFT");
//qa.checkForDefect("BSAWrong");
//qa.checkForDefect("BSAUnknown");
//qa.checkForDefect("TSAWrong");
//qa.checkForDefect("TSAUnknown");
//qa.checkForDefect("DSAWrong");
//qa.checkForDefect("DSAUnknown");
qa.checkForDefect("ChargeHigh");
qa.checkForDefect("ChargeNegative");
qa.checkForDefect("ChargeUnknown");
qa.checkForDefect("PossiblyNoBeam");

// ---------- iterate only over allowed runs ----------
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

// ---------- totals ----------
println "\ntotal accumulated charge: " + qa.getAccumulatedCharge() + " (nC)"
println "total no hel accumulated charge: " + qa.getAccumulatedChargeHL(0) + " (nC)"
println "total pos hel accumulated charge: " + qa.getAccumulatedChargeHL(1) + " (nC)"
println "total neg hel accumulated charge: " + qa.getAccumulatedChargeHL(-1) + " (nC)"
println "total unassigned hel accumulated charge: " + qa.getAccumulatedChargeHL(0) + " (nC)"
