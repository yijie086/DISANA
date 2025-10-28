// calculate total analyzed charge for a run period,
// with specified QA cuts enabled
// note: if syncCheck.groovy errors are present in the run range,
// the final charge value might be a bit wrong...

import org.jlab.io.hipo.HipoDataSource
import clasqa.QADB


// arguments: run range, from runLB to runUB
def runLB,runUB
if(args.length==2) {
  runLB = args[0].toInteger()
  runUB = args[1].toInteger()
}
else { System.err << "ARGUMENTS: [runLB] [runUB]\n"; return; }
println "run range: $runLB to $runUB"
QADB qa = new QADB("latest",runLB,runUB)
qa.checkForDefect("TotalOutlier");
qa.checkForDefect("TerminalOutlier");
qa.checkForDefect("MarginalOutlier");
qa.checkForDefect("SectorLoss");
qa.checkForDefect("LowLiveTime");
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


// loop over runs and files
int runnum
int evnum
qa.getQaTree().each{ runnumStr, runTree ->
  runnum = runnumStr.toInteger()
  // snapshot before this run
  double beforeTotal = qa.getAccumulatedCharge()
  double beforePos   = qa.getAccumulatedChargeHL(1)
  double beforeNeg   = qa.getAccumulatedChargeHL(-1)
  double before0   = qa.getAccumulatedChargeHL(0)
  runTree.each{ filenumStr, fileTree ->
    evnum = fileTree['evnumMin']

    // QA cut /////////////////
    //if(qa.query(runnum,evnum)) qa.accumulateCharge() // no qa cut
    //if(qa.golden(runnum,evnum)) qa.accumulateCharge()
    //if(qa.OkForAsymmetry(runnum,evnum)) qa.accumulateCharge()
    if(qa.pass(runnum, evnum)) {
    	qa.accumulateCharge();
    	qa.accumulateChargeHL();
    }
  }
  // calculate this run's charge
  double runTotal = qa.getAccumulatedCharge() - beforeTotal
  double runPos   = qa.getAccumulatedChargeHL(1) - beforePos
  double runNeg   = qa.getAccumulatedChargeHL(-1) - beforeNeg
  double run0   = qa.getAccumulatedChargeHL(0) - before0

  // print per‐run
  println "${runnum},${runTotal},${runPos},${runNeg},${run0}"
}

// print charge
println "\ntotal accumulated charge: " + qa.getAccumulatedCharge() + " (nC)"
println "total no hel accumulated charge: " + qa.getAccumulatedChargeHL(0) + " (nC)"
println "total pos hel accumulated charge: " + qa.getAccumulatedChargeHL(1) + " (nC)"
println "total neg hel accumulated charge: " + qa.getAccumulatedChargeHL(-1) + " (nC)"
println "total unassigned hel accumulated charge: " + qa.getAccumulatedChargeHL(0) + " (nC)"

