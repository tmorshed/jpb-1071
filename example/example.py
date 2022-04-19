from neuron import h,gui
from random import random as rand
from random import choice as pick

h.load_file("./GPe_Template.hoc")
h.load_file("./STN_Template.hoc")

numSTN = 1
numGPe = 1

STNs = []
GPes = []

for i in range(numSTN):
    STNs.append(h.STN())
for i in range(numGPe):
    GPes.append(h.GPe())

iclamp       = h.IClamp(STNs[0].soma[0](0.5))
iclamp.dur   = 0.1 # ms
iclamp.amp   = 15  # nA
iclamp.delay = 50  # ms

# STN Synapses [Inhibitory]

STN_syns = []
for stn in STNs:
    syn = h.ExpSyn(pick(stn.dend)(0.5))
    syn.e    = -80 # Inhibitory
    syn.tau  = 30 # Decay
    STN_syns.append(syn)

# GPe Synapses [Excitatory]

GPe_syns = []
for gpe in GPes:
    syn = h.ExpSyn(pick(gpe.dend)(0.5))
    syn.e    = 0 # Excitatory
    syn.tau  = 10 # Decay
    GPe_syns.append(syn)

# Connectivity

netcons = []
    
# STN to GPe Connections

for stn in STNs:
    axon = pick(stn.node)
    con = h.NetCon(axon(0.5)._ref_v,pick(GPe_syns),sec=axon)
    con.weight[0] = 10
    netcons.append(con)

# GPe to STN Connections

for gpe in GPes:
    axon = pick(gpe.node)
    con = h.NetCon(axon(0.5)._ref_v,pick(STN_syns),sec=axon)
    con.weight[0] = 10
    netcons.append(con)

# Recording (to be done before simulation)

STN_somas = []
GPe_somas = []

for stn in STNs:
    STN_somas.append(h.Vector().record(stn.soma[0](0.5)._ref_v))
for gpe in GPes:
    GPe_somas.append(h.Vector().record(gpe.soma[0](0.5)._ref_v))
    
# Simulation Control

h.finitialize(-65) # Initialize the simulator and
                   # all neuron sections to -65mV

tstop = 100 # ms
                   
while h.t < tstop:
    h.fadvance()
    if (int(h.t/tstop*1000)%10==0):
        print(f"{h.t/tstop*100:.1f}%\r",end="")
print("Done!")

import matplotlib.pyplot as plt

plt.subplot(121)
for stn in STN_somas:
    plt.plot(stn)
plt.subplot(122)
for gpe in GPe_somas:
    plt.plot(gpe)
plt.show()
