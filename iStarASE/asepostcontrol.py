rf = (1-918.0/1028.0)

def  asepostcontrol(thck,topg,Cwshelf, C, muCoef, divuho, x, y, *etc):

    #reground the deep hole behind smith
    if ( (x > 330.0e+3) & (x < 360.0e+3) & (y > 220e+3) &  (y < 245e+3)):
        usrf = max(thck + topg, rf*thck)
        Hf = usrf / rf
        topg = max(topg, Hf * (rf - 1.0) )
        divuho = 0.0

    return thck,topg,Cwshelf, C, muCoef, divuho
