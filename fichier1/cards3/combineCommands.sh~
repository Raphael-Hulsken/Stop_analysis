rm allLimits.root
rm allLimits2.root
rm higgsCombine*.root
combine -M ProfileLikelihood --significance 300_0.tab -t -1 --expectSignal=1 --toysFreq -n 300_0
combine -M ProfileLikelihood --significance 300_50.tab -t -1 --expectSignal=1 --toysFreq -n 300_50
combine -M ProfileLikelihood --significance 300_100.tab -t -1 --expectSignal=1 --toysFreq -n 300_100
combine -M ProfileLikelihood --significance 350_0.tab -t -1 --expectSignal=1 --toysFreq -n 350_0
combine -M ProfileLikelihood --significance 350_50.tab -t -1 --expectSignal=1 --toysFreq -n 350_50
combine -M ProfileLikelihood --significance 350_100.tab -t -1 --expectSignal=1 --toysFreq -n 350_100
combine -M ProfileLikelihood --significance 350_150.tab -t -1 --expectSignal=1 --toysFreq -n 350_150
combine -M ProfileLikelihood --significance 400_0.tab -t -1 --expectSignal=1 --toysFreq -n 400_0
combine -M ProfileLikelihood --significance 400_25.tab -t -1 --expectSignal=1 --toysFreq -n 400_25
combine -M ProfileLikelihood --significance 400_50.tab -t -1 --expectSignal=1 --toysFreq -n 400_50
combine -M ProfileLikelihood --significance 400_100.tab -t -1 --expectSignal=1 --toysFreq -n 400_100
hadd allLimits.root higgsCombine*.Asymptotic.mH120.root
hadd allLimits2.root higgsCombine*.ProfileLikelihood.mH120.root
