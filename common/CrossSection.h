#ifndef CrossSection_h
#define CrossSection_h


float CrossSection(const string& dataset){
	if(dataset == "T2tt_425_325") return 1.31169;
	if(dataset == "T2tt_500_325") return 0.51848;
	if(dataset == "T2tt_650_100") return 0.107045;
	if(dataset == "T2tt_650_325") return 0.107045;
	if(dataset == "T2tt_850_100") return 0.01896;
        if(dataset == "T2tt_850_325") return 0.01896;
	if(dataset == "singleTopbar_s") return 4.16;
	if(dataset == "singleTopbar_t") return 80.95;
	if(dataset == "singleTop_s") return 7.20;
	if(dataset == "singleTop_t") return 136.02;
	if(dataset == "ttbar-madgraph") return 831.76;
	if(dataset == "TTW") return 0.70;
	if(dataset == "TTW_ln") return 0.70*0.32;
	if(dataset == "TTW_qq") return 0.70*0.675;
	if(dataset == "TTZ") return 0.62;
	if(dataset == "TTZ_ll") return 0.62*0.2;
	if(dataset == "TTZ_qq") return 0.62*0.7;
	if(dataset == "Wjets") return 61466;
	if(dataset == "WZ") return 48.4;
	if(dataset == "WW_aMC") return 48.4;
	if(dataset == "ZZ_aMC") return 15.4;
	if(dataset == "ZZ") return 15.4;
	if(dataset == "TTJets") return  831.76 ;
	if(dataset == "TTjets_M5") return  831.76 ;
        if(dataset == "TTJets-FXFX") return  831.76 ; 
        if(dataset == "WJets") return 61466;
        if(dataset == "Wjets_aMC") return 61466;
    	if(dataset == "ST_s") return  	10.32;
    	if(dataset == "ST_t-atop") return 80.95;
    	if(dataset == "ST_t-top") return  136.02 ;
    	if(dataset == "ST_tW-atop") return 35.6 ;
    	if(dataset == "ST_tW-top") return 35.6 ;
    return 0;
}

#endif
