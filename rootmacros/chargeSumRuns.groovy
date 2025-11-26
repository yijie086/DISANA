#!/usr/bin/env groovy
// Calculate total analyzed charge with QA cuts using only run numbers.
// Modes:
//   1) run-only:      runLB runUB
//   2) file-only:     [runListFile]
// Optional flags anywhere: --exclude "3211,3215-3220"  --exclude-file file
//
// This script assumes the QA DB is available for the specified runs.
// It removes all file system scanning (paths to .hipo files) and the --recursive flag.

import clasqa.QADB
import java.nio.file.*

// ---------- helpers ----------

boolean isInt(String s) {
  try { s?.toInteger(); return true } catch (ignored) { return false }
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

// Loads run numbers from a text file (like charge_sp19.txt)
Set<Integer> loadRunListFile(String path) {
  def out = new LinkedHashSet<Integer>()
  if (!path) return out
  def f = new File(path)
  if (!f.exists() || !f.isFile()) die("ERROR: run list file '${path}' not found.")
  println "Parsing runs from: ${path}"
  // This logic assumes the run number is the first token on a line and is 4+ digits.
  f.eachLine { line ->
    // This regex looks for 4 or more digits at the start of a line, optionally preceded by whitespace.
    def match = (line =~ /^\s*(\d{4,})/)
    if (match.find()) {
      try { out << match.group(1).toInteger() } catch (ignored) {}
    }
  }
  if (out.isEmpty()) die("ERROR: Could not parse any run numbers from '${path}'.")
  return out
}

// ---------- parse flags (can appear anywhere) ----------
def argList = (args as List).collect { it as String }
// Removed: boolean recursive = false

String excludeSpec = null
String excludeFile = null

for (int i = 0; i < argList.size();) {
  def a = argList[i]
  // Removed: if (a == "--recursive") { ... }
  if (a == "--exclude" && i+1 < argList.size()) { excludeSpec = argList[i+1]; argList.remove(i+1); argList.remove(i) }
  else if (a == "--exclude-file" && i+1 < argList.size()) { excludeFile = argList[i+1]; argList.remove(i+1); argList.remove(i) }
  else { i++ }
}

// ---------- decide mode from remaining positionals ----------
Integer runLB = null, runUB = null
File runListFile = null
boolean runOnlyMode = false
boolean fileOnlyMode = false

if (argList.size() == 0) {
  die("ERROR: No arguments provided. Use: [runLB runUB] or [runListFile]")
} else if (argList.size() == 1) {
  // single positional: a run list file
  if (isInt(argList[0])) {
    die("ERROR: Received one integer ('${argList[0]}'). For run-only mode provide: runLB runUB")
  }
  
  runListFile = new File(argList[0])
  fileOnlyMode = true 
  
} else if (argList.size() == 2 && isInt(argList[0]) && isInt(argList[1])) {
  // have runLB runUB
  runLB = argList[0].toInteger()
  runUB = argList[1].toInteger()
  if (runLB > runUB) die("ERROR: runLB (${runLB}) must be <= runUB (${runUB}).")
  runOnlyMode = true
} else {
  die("ERROR: Unrecognized arguments. Use either [runLB runUB] or [runListFile].")
}

// ---------- build allowedRuns ----------
Set<Integer> allowedRuns
if (runOnlyMode) {
  allowedRuns = (runLB..runUB).toSet()
  println "Mode: run-only. Using runs ${runLB}..${runUB}."
} else if (fileOnlyMode) {
  allowedRuns = loadRunListFile(runListFile.absolutePath)
  runLB = allowedRuns.min(); runUB = allowedRuns.max()
  println "Mode: file-only. Inferred run range from list file: ${runLB} to ${runUB}"
} 

// ---------- apply exclusions (all modes) ----------
def excludeRuns = new LinkedHashSet<Integer>()
excludeRuns.addAll(parseRunList(excludeSpec))
excludeRuns.addAll(loadExcludeFile(excludeFile))
excludeRuns.addAll(parseRunList(System.getenv("EXCLUDE_RUNS")))
if (!allowedRuns.isEmpty() && !excludeRuns.isEmpty()) {
  def before = allowedRuns.size()
  allowedRuns.removeAll(excludeRuns)
  println "Excluded ${excludeRuns.size()} specified runs; kept ${allowedRuns.size()} of ${before}."
  if (allowedRuns.isEmpty()) die("No runs left to analyze after exclusions.")
}

if (allowedRuns.isEmpty()) die("No runs found to analyze.")

println "Run range considered: $runLB to $runUB"
println "Number of runs to analyze: ${allowedRuns.size()}"

// ---------- QADB Configuration: Allow 'Misc' bit for specific runs ----------
// Initialize a temporary QA object covering the full range to apply configuration updates
//println "Configured QADB: Allowed 'Misc' bit for ${configQa.getAllowedMiscRuns().size()} specific runs."
// configQa instance is no longer needed after configuration is applied.

// ---------- QA setup (per-run helper) ----------
QADB makeQaForRun(int run) {
  QADB qa = new QADB("latest", run, run)

  // Enable desired defect checks
qa.checkForDefect('TotalOutlier')    
qa.checkForDefect('TerminalOutlier')
qa.checkForDefect('MarginalOutlier')
qa.checkForDefect('SectorLoss')
qa.checkForDefect('Misc')
qa.checkForDefect('ChargeHigh')
qa.checkForDefect('ChargeNegative')
qa.checkForDefect('ChargeUnknown')
qa.checkForDefect('PossiblyNoBeam')
  return qa
}

// ---------- iterate only over allowed runs (per-run QADB, manual totals) ----------
double totalAll       = 0.0
double totalAllPos    = 0.0
double totalAllNeg    = 0.0
double totalAllNoHel  = 0.0
double totalAllUnkHel = 0.0   // here we treat 0 helicity as unassigned

allowedRuns.each { int runnum ->

  QADB qa = makeQaForRun(runnum)
[ // List of runs with `Misc` that should be allowed
    5046, 5047, 5051, 5128, 5129, 5130, 5158, 5159,
    5160, 5163, 5165, 5166, 5167, 5168, 5169, 5180,
    5181, 5182, 5183, 5400, 5448, 5495, 5496, 5505,
    5567, 5610, 5617, 5621, 5623, 6736, 6737, 6738,
    6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747,
    6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756,
    6757,                        // RGA runs ^
      16194, 16089, 16185, 16308, 16184, 16307, 16309, // RGC Su22 He/ET
      16872, 16975,                    // RGC Fa22 He/ET
      17763, 17764, 17765, 17766, 17767, 17768,    // RGC Sp23 He/ET
    17179, 17180, 17181, 17182, 17183, 17188, 17189, // RICH off/partially down
      17252
].each{ run -> qa.allowMiscBit(run) }

  def qaTree = qa.getQaTree()
  if (qaTree == null || qaTree.isEmpty()) {
    System.err.println("WARNING: no QA tree for run ${runnum}")
    return
  }

  qaTree.each { runnumStr, runTree ->
    int rn = runnumStr.toInteger()
    runTree.each { filenumStr, fileTree ->
      int evnum = fileTree['evnumMin']
      if (qa.pass(rn, evnum)) {
        qa.accumulateCharge()
        qa.accumulateChargeHL()
      }
    }
  }

  double runTotal = qa.getAccumulatedCharge()
  double runPos   = qa.getAccumulatedChargeHL(1)
  double runNeg   = qa.getAccumulatedChargeHL(-1)
  double run0     = qa.getAccumulatedChargeHL(0)

  println "${runnum},${runTotal},${runPos},${runNeg},${run0}"

  totalAll       += runTotal
  totalAllPos    += runPos
  totalAllNeg    += runNeg
  totalAllNoHel  += run0
  totalAllUnkHel += run0
}

// ---------- totals ----------
println "\ntotal accumulated charge: " + totalAll + " (nC)"
println "total no hel accumulated charge: " + totalAllNoHel + " (nC)"
println "total pos hel accumulated charge: " + totalAllPos + " (nC)"
println "total neg hel accumulated charge: " + totalAllNeg + " (nC)"
println "total unassigned hel accumulated charge: " + totalAllUnkHel + " (nC)"