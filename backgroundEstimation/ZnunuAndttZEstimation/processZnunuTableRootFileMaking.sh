make ZnunuUnctTablePerProcess
make ZnunuUnctRootFilePerProcess
make ZnunuUnctSummedTables
./ZnunuUnctTablePerProcess yieldMorttZ.tab signalRegMor.txt statNames.txt 17 yieldMorZNuNuJECDown.tab yieldMorZNuNuJECUp.tab ttZ notgroup
./ZnunuUnctTablePerProcess yieldMorWZ.tab signalRegMor.txt statNames.txt 17 yieldMorZNuNuJECDown.tab yieldMorZNuNuJECUp.tab WZ notgroup
./ZnunuUnctTablePerProcess yieldMorZZ.tab signalRegMor.txt statNames.txt 17 yieldMorZNuNuJECDown.tab yieldMorZNuNuJECUp.tab ZZ notgroup
./ZnunuUnctRootFilePerProcess ttZtableUnc.tab realregions.txt uncertainties.txt ttZ
./ZnunuUnctRootFilePerProcess WZtableUnc.tab realregions.txt uncertainties.txt WZ
./ZnunuUnctRootFilePerProcess ZZtableUnc.tab realregions.txt uncertainties.txt ZZ
./ZnunuUnctSummedTables ttZtableUnc.tab WZtableUnc.tab ZZtableUnc.tab realregions.txt uncertainties.txt notgroup

