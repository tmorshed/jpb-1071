TITLE calcium accumulation for GPi neuron model

COMMENT 

    From: Gillies and Willshaw (2006) and as implemented by Miocinovic (2006)

    Calcium accumulation into a volume of area*depth next to the
    membrane with an exponential decay (time constant tau) to resting
    level (given by the global calcium variable cai0_ca_ion).

ENDCOMMENT

NEURON {
    SUFFIX GPi_cacum
    USEION ca READ ica WRITE cai, cao
    GLOBAL depth, tau, cai0, cao0
}

UNITS {
    (mM) = (milli/liter)
    (mA) = (milliamp)
    F = (faraday) (coulombs)	: Faraday constant 
}

PARAMETER {
    depth = 200 (nm)      : assume volume = area * depth
    tau = 10 (ms)
    cai0 = 5e-5 (mM)
	cao0 = 2.4	(mM)
    cai0_ca_ion = 1e-4 (mM)
}

ASSIGNED {
    ica (mA/cm2)
}

STATE {
    cai (mM)
	cao	(mM)
}

INITIAL {
    cai=cai0
	cao=cao0
}

BREAKPOINT {
    SOLVE integrate METHOD derivimplicit
}

DERIVATIVE integrate {
    cai' = -ica/depth/F/2 * (1e7) + (cai0 - cai)/tau
	cao' = 0
}
