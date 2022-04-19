#%% Import and load
from neuron import h,gui
from random import random as rand
from random import choice as pick
from random import seed as seed
from tqdm.notebook import tqdm as tqdm
import numpy as np
import matplotlib.pyplot as plt

h.load_file("./GPe_Template.hoc")
h.load_file("./STN_Template.hoc")

# set random seed
# seed(10)

#%% Inititate

def simulate(stn_con_wt=10, gpe_con_wt = 15, numSTN =3, numGPe = 3, tstop = 100):
    #tstop: ms
    STNs = []
    GPes = []

    # for i in range(numSTN):
    #     STNs.append(h.STN())
    # for i in range(numGPe):
    #     GPes.append(h.GPe())

    STNs = [h.STN() for i in range(numSTN)]
    GPes = [h.GPe() for i in range(numGPe)]

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
        con.weight[0] = stn_con_wt
        netcons.append(con)

    # GPe to STN Connections

    for gpe in GPes:
        axon = pick(gpe.node)
        con = h.NetCon(axon(0.5)._ref_v,pick(STN_syns),sec=axon)
        con.weight[0] = gpe_con_wt
        netcons.append(con)

    # Recording (to be done before simulation)

    STN_somas = []
    GPe_somas = []

    for stn in STNs:
        STN_somas.append(h.Vector().record(stn.soma[0](0.5)._ref_v))
    for gpe in GPes:
        GPe_somas.append(h.Vector().record(gpe.soma[0](0.5)._ref_v))

    # simulation:
    h.finitialize(-65) # Initialize the simulator and
                    # all neuron sections to -65mV


    for i in tqdm(np.arange(h.t,tstop,h.dt,dtype=float), desc="single sim"):
        h.fadvance()
    return STN_somas, GPe_somas, numSTN, numGPe, tstop

# while h.t < tstop:
#     h.fadvance()
#     if (int(h.t/tstop*1000)%10==0):
#         print(f"{h.t/tstop*100:.1f}%\r",end="")
# print("Done!")

#%% Visualize

def visualize(STN_somas, GPe_somas):
    nrows=max([numSTN, numGPe])
    fig, ax = plt.subplots(figsize=(15,7.5), ncols=2, nrows=nrows, sharey=True)
    
    # STNs visualization:
    for i,j in enumerate(ax[:,0]):
        j.plot(np.linspace(0,tstop, len(STN_somas)), STN_somas[:,i])
        j.set(title="STN ["+str(i)+"]", ylabel="Voltage (mV)", xlabel='time (ms)')
    
    for i,j in enumerate(ax[:,1]):
        j.plot(np.linspace(0,tstop, len(GPe_somas)), GPe_somas[:,i])
        j.set(title="GPe ["+str(i)+"]", ylabel="Voltage (mV)", xlabel='time (ms)')

    # ax[1].plot(np.linspace(0,tstop, len(GPe_somas)), GPe_somas)
    # ax[1].set(title="GPe", ylabel="Voltage (mV)", xlabel='time (ms)')

    plt.suptitle("Connectivity weights: STN =" + str(stn_con_wt) + " & GPe =" + str(gpe_con_wt))
    plt.tight_layout()
    plt.savefig("./figs/" + str(stn_con_wt) + "s" + str(gpe_con_wt) + "g")
# %% RUN SINGLE SIM

stn_con_wt = 10
gpe_con_wt = 15
numSTN = 3
numGPe = 3
tstop = 100

STN_somas=np.array(simvals[0]).T
GPe_somas=np.array(simvals[1]).T

simvals = simulate(stn_con_wt=stn_con_wt,gpe_con_wt=gpe_con_wt, numSTN=numSTN, numGPe=numGPe)
visualize(STN_somas,GPe_somas)
# %% RUN MULTIPLE SIMS

stn_con_wts = [10,15,20]
gpe_con_wts = [10,15,20]
numSTN = 3
numGPe = 3
tstop = 100

with tqdm(total=9,colour='blue') as pbar:
    pbar.set_description("TOTAL Progress")
    for stn_con_wt in stn_con_wts:
        for gpe_con_wt in gpe_con_wts:
            simvals = simulate(stn_con_wt=stn_con_wt,gpe_con_wt=gpe_con_wt, numSTN=numSTN, numGPe=numGPe)
            STN_somas=np.array(simvals[0]).T
            GPe_somas=np.array(simvals[1]).T
            visualize(STN_somas,GPe_somas)
            pbar.update(1)
# %%
