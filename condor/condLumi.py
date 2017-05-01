#MuRunBv3
crabLumiBv3 = 5.19
condLumiBv3 =  (60.0/62.0)
finalLumiBv3 = crabLumiBv3*condLumiBv3
print "% of condor success: ", condLumiBv3*100
print "Bv3 lumi after cond: ", finalLumiBv3


#MuRunCv1
crabLumiCv1 = 2.25
condLumiCv1 =  (33.0/35)
finalLumiCv1 = crabLumiCv1*condLumiCv1
print "\n % of condor success: ", condLumiCv1*100
print "Cv1 lumi after cond: : ", finalLumiCv1

#MuRunDv1
crabLumiDv1 = 4.25
condLumiDv1 =  (51.0/53)
finalLumiDv1 = crabLumiDv1*condLumiDv1
print "\n % of condor success: ", condLumiDv1*100
print "Dv1 lumi after cond: : ", finalLumiDv1

#MuRunEv1
crabLumiEv1 = 3.95
condLumiEv1 =  (42.0/44)
finalLumiEv1 = crabLumiEv1*condLumiEv1
print "\n % of condor success: ", condLumiEv1*100
print "Ev1 lumi after cond: : ", finalLumiEv1


#MuRunFv1
crabLumiFv1 = 2.94
condLumiFv1 =  (36.0/37)
finalLumiFv1 = crabLumiFv1*condLumiFv1
print "\n % of condor success: ", condLumiFv1*100
print "Fv1 lumi after cond: : ", finalLumiFv1


#MuRunGv1
crabLumiGv1 = 7.43
condLumiGv1 =  (43.0/71)
finalLumiGv1 = crabLumiGv1*condLumiGv1
print "\n % of condor success: ", condLumiGv1*100
print "Gv1 lumi after cond: : ", finalLumiGv1


#MuRunHv2
crabLumiHv2 = 8.15
condLumiHv2 =  (81.0/88)
finalLumiHv2 = crabLumiHv2*condLumiHv2
print "\n % of condor success: ", condLumiHv2*100
print "Hv2 lumi after cond: : ", finalLumiHv2


#MuRunHv3
crabLumiHv3 = 0.22
condLumiHv3 =  (03.0/03)
finalLumiHv3 = crabLumiHv3*condLumiHv3
print "\n % of condor success: ", condLumiHv3*100
print "Hv3 lumi after cond: : ", finalLumiHv3

totCrabLumi = crabLumiBv3 + crabLumiCv1 + crabLumiDv1+ crabLumiEv1+ crabLumiFv1+ crabLumiGv1+ crabLumiHv2+ crabLumiHv3
totCondLumi = finalLumiBv3 + finalLumiCv1 + finalLumiDv1 + finalLumiEv1 + finalLumiFv1 + finalLumiGv1 + finalLumiHv2+ finalLumiHv3

print "-------------------------------------"
print "Total Crab Luminosity: ", totCrabLumi
print "Total Cond luminosity: ", totCondLumi
print "-------------------------------------"








