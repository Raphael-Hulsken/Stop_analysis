#ifndef SIGNAL_CS
#define SIGNAL_CS
#include <exception>
float returnSigCS(float stopm)
{
    if(stopm == 100 ) return 1521.11 ;
    else if(stopm == 125) return 574.981;
    else if(stopm == 150) return 249.409;
    else if(stopm == 175) return 121.416;
    else if(stopm == 200) return 64.5085;
    else if(stopm == 225) return 36.3818;
    else if(stopm == 250) return 21.5949;
    else if(stopm == 275) return 13.3231;
    else if(stopm == 300) return 8.51615;
    else if(stopm == 325) return 5.60471;
    else if(stopm == 350) return 3.78661;
    else if(stopm == 375) return 2.61162;
    else if(stopm == 400) return 1.83537;
    else if(stopm == 425) return 1.31169;
    else if(stopm == 450) return 0.948333;
    else if(stopm == 475) return 0.697075;
    else if(stopm == 500) return 0.51848;
    else if(stopm == 525) return 0.390303;
    else if(stopm == 550) return 0.296128;
    else if(stopm == 600) return 0.174599;
    else if(stopm == 650) return 0.107045;
    else if(stopm == 700) return 0.0670476;
    else if(stopm == 750) return 0.0431418;
    else if(stopm == 800) return 0.0283338;
    else if(stopm == 850) return 0.0189612;
    else if(stopm == 900) return 0.0128895;
    else if(stopm == 950) return 0.00883465;
    else if(stopm == 1000) return 0.00615134;
    else if(stopm == 1050) return 0.00432261;
    else if(stopm == 1100) return 0.00307413;
    else if(stopm == 1150) return 0.00221047;
    else if(stopm == 1200) return 0.00159844;
    else if(stopm == 1250) return 0.0011583;
    else 
    {
        std::cout << "stop mass: " << stopm << std::endl;
        std::cout << "cheating a lot, the stop xsec not correct, this is just for testing purposes" << endl;
        //throw std::runtime_error("no cross section for given stop mass"); //@MJ@ TODO make this working!
        return 0.1;
    }

}

#endif
