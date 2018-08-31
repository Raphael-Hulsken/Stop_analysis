rm allLimits.root
rm allLimits2.root
rm higgsCombine*.root
combine -M Asymptotic 500_300.tab -n 500_300
combine -M ProfileLikelihood --significance 500_300.tab -t -1 --expectSignal=1 --toysFreq -n 500_300
hadd allLimits.root higgsCombine*.Asymptotic.mH120.root
hadd allLimits2.root higgsCombine*.ProfileLikelihood.mH120.root
