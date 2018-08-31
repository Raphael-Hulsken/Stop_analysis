rm allLimits.root
rm allLimits2.root
rm higgsCombine*.root
combine -M Asymptotic 625_425.tab -n 625_425
combine -M ProfileLikelihood --significance 625_425.tab -t -1 --expectSignal=1 --toysFreq -n 625_425
hadd allLimits.root higgsCombine*.Asymptotic.mH120.root
hadd allLimits2.root higgsCombine*.ProfileLikelihood.mH120.root
