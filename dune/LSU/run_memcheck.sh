#/bin/bash

valgrind --track-origins=yes --read-var-info=yes --partial-loads-ok=no --num-callers=25 `which lar` -n 5 -c piAbsSelector_redo_beamevent.fcl /cshare/vol1/data/protodune/reco_v07_08_04/run5145-7GeV-180kV/np04_raw_run005145_0001_dl1_reco_13141411_0_20181110T060526.root -T dumb.root >& memcheck.txt
