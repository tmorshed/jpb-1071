TITLE default extra/intra-cellular ion concentrations for GPe

COMMENT
 Used to set: nai, nao, ki, ko
   so that NEURON's default values are not used

ENDCOMMENT

NEURON {
    SUFFIX GPi_myions
    USEION k WRITE ko, ki
    USEION na WRITE nao, nai
	RANGE nai0, nao0, ki0, ko0
}

UNITS {
    (molar)	= (1/liter)
    (mM)	= (millimolar)
}

PARAMETER {
	nai0 = 10	(mM)
	nao0 = 150	(mM)
	ki0  = 140	(mM)
	ko0  = 6.24	(mM)
}
	
STATE {
    nai (mM)
    nao (mM)
    ki (mM)
    ko (mM)
}

INITIAL {
    nai = nai0
    nao = nao0
    ki = ki0
    ko = ko0
}