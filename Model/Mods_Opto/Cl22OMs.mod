COMMENT

Cl22OMs: chloride conducting 22OMs opsin model

ENDCOMMENT

NEURON{
    SUFFIX Cl22OMs
    POINTER Iopto
    USEION cl READ cli, clo, ecl WRITE icl VALENCE -1
    RANGE gmax, g, o, da, Erev, tauo, tauda, oinf, dainf, Rinf, EvarFlag, iopsin
    RANGE Iratio, tauOn, tauOff, tauInact, tauRecov, Oinf_p1, Oinf_p2
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mS) = (millisiemens)
  (um) = (micron)
  (S)  = (siemens)
  (W)  = (watt)
} 

PARAMETER{

    gmax = 0.012885323820985 (S/cm2)
    Erev = -85 (mV)

    EvarFlag = 0

    clo
    cli

    dainf
    oinf

    scaling

    Oinf_p1 = 2.837 : 2.473
    Oinf_p2 = 0.942 : 0.677

    Iratio = 0.8

    tauOn = 0.018
    tauOff = 140.014
    tauInact = 236.82
    tauRecov = 4803
    
}

ASSIGNED{
    v (mV)
    icl (mA/cm2)
    ecl (mV)
    g (S/cm2)
    tauo (s)
    tauda (s)

    Iopto (W/m2)

    iopsin (mA/cm2)
    
}

STATE{
    o da
}

INITIAL{    
    o=0
    da=1
    
    Iopto=0

    scaling = 0
    
}

BREAKPOINT{
    
    SOLVE states METHOD cnexp

    g = gmax*o*da
    if (EvarFlag == 0) {iopsin = g*(v-Erev)}
    if (EvarFlag == 1) {iopsin = g*(v-ecl)}
    
    icl = iopsin
}

DERIVATIVE states{
    ss(-65, Iopto)
    o' = (oinf-o)/tauo
    da' = (dainf-da)/tauda
}

PROCEDURE ss(vm(mV),Io(W/m2)){
    if (Io == 0) {oinf = 0} else {oinf = 1/(1+exp(Oinf_p1/Oinf_p2)*Iopto^(-1/(Oinf_p2*log(10))))}
    if (Io == 0) {dainf = 1} else {dainf = Iratio}
    if (Io == 0) {tauo = tauOff} else {tauo = tauOn}
    if (Io == 0) {tauda = tauRecov} else {tauda = tauInact}
}
