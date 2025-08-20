COMMENT

K22OM: potassium conducting 22OM opsin model

ENDCOMMENT

NEURON{
    SUFFIX K22OM
    POINTER Iopto
    USEION k READ ki, ko, ek WRITE  ik VALENCE 1
    RANGE gmax, g, o, da, Erev, Iintm, tauo, tauda, oinf, dainf, EvarFlag, iopsin
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mS) = (millisiemens)
  (um) = (micron)
  (S)     = (siemens)
  (W) = (watt)
} 

PARAMETER{

    gmax = 0.012885323820985 (S/cm2)
    Erev = -69.5 (mV)

    EvarFlag = 0

    ko
    ki

    Oinf_p1 = 2.473
    Oinf_p2 = 0.677

    Rinf_p1 = 2.364
    Rinf_p2 = 0.195
    Rinf_p3 = 0.323

    tauOi_p1 = 0.6
    tauOi_p2 = -1.007
    tauOi_p3 = 0.088

    tauRi_p1 = 0.194
    tauRi_p2 = 0.047
    tauRi_p3 = 2.453
    tauRi_p4 = 0.154
    tauRi_p5 = -0.699
    tauRi_p6 =  1.522

    tauOv_p1 = 0.879
    tauOv_p2 = 80.064
    tauOv_p3 = -72.976

    tauRv_p1 = 63.095
    tauRv_p2 = -96.234
    tauRv_p3 = -95.040
}

ASSIGNED{
    v (mV)
    ik (mA/cm2)
    ek (mV)
    vk (mV)
    g (S/cm2)
    tauo (s)
    tauda (s)
	oinf
	dainf

    Iopto (W/m2)
    Iintm (W/m2)

    iopsin (mA/cm2)
}

STATE{
    o da
}

INITIAL{    
    o=0
    da=1
    
    Iopto=0
    
}

BREAKPOINT{
    SOLVE states METHOD cnexp
    g = gmax*o*da
    if (EvarFlag == 0) {iopsin = g*(v-Erev)}
    if (EvarFlag == 1) {iopsin = g*(v-ek)}
    
    ik = iopsin
    vk = v
}

DERIVATIVE states{
    ss(v,Iopto)   

    o' = (oinf-o)/tauo
    da' = (dainf-da)/tauda
    
}


PROCEDURE ss(vm(mV),Io(W/m2)){

    Iintm = Io

    oinf = 1/(1+exp(Oinf_p1/Oinf_p2)*Io^(-1/(Oinf_p2*log(10))))

    dainf = 1-(Rinf_p3/(1+exp(Rinf_p1/Rinf_p2)*Io^(-1/(Rinf_p2*log(10)))))
    
    tauo =tauOv_p1/(1+exp(-(vm-tauOv_p2)/tauOv_p3))*tauOi_p3/(1+exp(tauOi_p1/tauOi_p2)*Io^(-1/(tauOi_p2*log(10))))*1e3

    tauda = tauRv_p1/(1+exp(-(vm-tauRv_p2)/tauRv_p3))*tauRi_p1*(1-tauRi_p2/(1+exp(tauRi_p3/tauRi_p4)*Io^(-1/(tauRi_p4*log(10))))-(1-tauRi_p2)/(1+exp(tauRi_p5/tauRi_p6)*Io^(-1/(tauRi_p6*log(10)))))*1e3

}

