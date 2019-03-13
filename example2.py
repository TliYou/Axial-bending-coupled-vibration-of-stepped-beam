from pylab import *
from beam import *


Es = 210e9
Ea = 69e9
rhos = 7800
rhoa = 2700
nus = 0.3
nua = 0.33
Ec = 110e9
nuc = 0.32
rhoc = 8500

#sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
#                [rhoc,rhoa], [0.85,0.85], 'section-1')
#sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
#                [rhoc,rhoa], [0.85,0.85], 'section-2')
#sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
#                [rhoc,rhoa], [0.85,0.85], 'section-3')
#beam = Beam()
#beam.append_segment(1.0, sect1)
#beam.append_segment(2.0, sect2)
#beam.append_segment(3.0, sect3)
#beam.nodes[0].add_spring(5e7,2e8)
#beam.nodes[1].add_mass(50)
#beam.nodes[2].add_mass(50)
#beam.nodes[3].add_spring(5e7,2e8)


def c_ku():
    ks = logspace(2,10,50)
    nfs = []
    for k in ks:
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam()
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].fixed()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].add_spring(k,0,0)
        f = beam.natural_frequency_real(0,500,200)[:5]
        nfs.append(f)
    nfs = np.array(nfs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    # adjust manually
#    j = 30
#    i,i2 = 1,3
#    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
#                       nfs[i,j+1:].copy())

    styles = ['-',':','--','.-','-.']
    labels = 'mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5'
    fig = figure(figsize=(4,3.6))
    ax = fig.add_subplot(111)
    for i in range(4):
        ax.loglog(ks,nfs[i],styles[i],label=labels[i])
    ax.set_xlabel(r'$k_{u2}\ (\mathrm{N/m})$', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
    ax.legend(loc='best')
#    xticks(fontsize=15)
#    yticks(fontsize=15)
    fig.tight_layout()
#    ax.set_ylim(0,300)
    savefig('example2_c_ku.png',dpi=800)
    savefig('example2_c_ku.svg',dpi=1000)
    return nfs

def c_kv():
    ks = logspace(2,10,50)
    nfs = []
    for k in ks:
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam()
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].fixed()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].add_spring(0,k,0)
        f = beam.natural_frequency_real(0,500,200)[:5]
        nfs.append(f)
    nfs = np.array(nfs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    # adjust manually
#    j = 30
#    i,i2 = 1,3
#    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
#                       nfs[i,j+1:].copy())

    styles = ['-',':','--','.-','-.']
    labels = 'mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5'
    fig = figure(figsize=(4,3.6))
    ax = fig.add_subplot(111)
    for i in range(4):
        ax.loglog(ks,nfs[i],styles[i],label=labels[i])
    ax.set_xlabel(r'$k_{v2}\ (\mathrm{N/m})$', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
    ax.legend(loc='best')
#    xticks(fontsize=15)
#    yticks(fontsize=15)
    fig.tight_layout()
#    ax.set_ylim(0,300)
    savefig('example2_c_kv.png',dpi=800)
    savefig('example2_c_kv.svg',dpi=1000)
    return nfs


def c_kr():
    ks = logspace(2,10,50)
    nfs = []
    for k in ks:
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam()
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].fixed()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].add_spring(0,0,k)
        f = beam.natural_frequency_real(0,500,200)[:5]
        nfs.append(f)
    nfs = np.array(nfs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    # adjust manually
#    j = 30
#    i,i2 = 1,3
#    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
#                       nfs[i,j+1:].copy())

    styles = ['-',':','--','.-','-.']
    labels = 'mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5'
    fig = figure(figsize=(4,3.6))
    ax = fig.add_subplot(111)
    for i in range(4):
        ax.loglog(ks,nfs[i],styles[i],label=labels[i])
    ax.set_xlabel(r'$k_{\theta 2}\ (\mathrm{N/m})$', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
    ax.legend(loc='best')
#    xticks(fontsize=15)
#    yticks(fontsize=15)
    fig.tight_layout()
#    ax.set_ylim(0,300)
    savefig('example2_c_kr.png',dpi=800)
    savefig('example2_c_kr.svg',dpi=1000)
    return nfs



def s_ku():
    ks = logspace(2,10,50)
    nfs = []
    for k in ks:
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam()
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].pinned()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].add_spring(k,0,0)
        f = beam.natural_frequency_real(0,500,200)[:5]
        nfs.append(f)
    nfs = np.array(nfs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    # adjust manually
#    j = 30
#    i,i2 = 1,3
#    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
#                       nfs[i,j+1:].copy())

    styles = ['-',':','--','.-','-.']
    labels = 'mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5'
    fig = figure(figsize=(4,3.6))
    ax = fig.add_subplot(111)
    for i in range(4):
        ax.loglog(ks,nfs[i],styles[i],label=labels[i])
    ax.set_xlabel(r'$k_{u2}\ (\mathrm{N/m})$', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
    ax.legend(loc='best')
#    xticks(fontsize=15)
#    yticks(fontsize=15)
    fig.tight_layout()
#    ax.set_ylim(0,300)
    savefig('example2_s_ku.png',dpi=800)
    savefig('example2_s_ku.svg',dpi=1000)
    return nfs

def s_kv():
    ks = logspace(2,10,50)
    nfs = []
    for k in ks:
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam()
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].pinned()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].add_spring(0,k,0)
        f = beam.natural_frequency_real(0,500,200)[:5]
        nfs.append(f)
    nfs = np.array(nfs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    # adjust manually
#    j = 30
#    i,i2 = 1,3
#    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
#                       nfs[i,j+1:].copy())

    styles = ['-',':','--','.-','-.']
    labels = 'mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5'
    fig = figure(figsize=(4,3.6))
    ax = fig.add_subplot(111)
    for i in range(4):
        ax.loglog(ks,nfs[i],styles[i],label=labels[i])
    ax.set_xlabel(r'$k_{v2}\ (\mathrm{N/m})$', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
    ax.legend(loc='best')
#    xticks(fontsize=15)
#    yticks(fontsize=15)
    fig.tight_layout()
#    ax.set_ylim(0,300)
    savefig('example2_s_kv.png',dpi=800)
    savefig('example2_s_kv.svg',dpi=1000)
    return nfs


def s_kr():
    ks = logspace(2,10,50)
    nfs = []
    for k in ks:
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam()
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].pinned()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].add_spring(0,0,k)
        f = beam.natural_frequency_real(0,500,200)[:5]
        nfs.append(f)
    nfs = np.array(nfs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    # adjust manually
#    j = 30
#    i,i2 = 1,3
#    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
#                       nfs[i,j+1:].copy())

    styles = ['-',':','--','.-','-.']
    labels = 'mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5'
    fig = figure(figsize=(4,3.6))
    ax = fig.add_subplot(111)
    for i in range(4):
        ax.loglog(ks,nfs[i],styles[i],label=labels[i])
    ax.set_xlabel(r'$k_{\theta 2}\ (\mathrm{N/m})$', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
    ax.legend(loc='best')
#    xticks(fontsize=15)
#    yticks(fontsize=15)
    fig.tight_layout()
#    ax.set_ylim(0,300)
    savefig('example2_s_kr.png',dpi=800)
    savefig('example2_s_kr.svg',dpi=1000)
    return nfs


def f_ku():
    ks = logspace(2,10,50)
    nfs = []
    for k in ks:
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam()
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].free()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].add_spring(k,0,0)
        f = beam.natural_frequency_real(0,500,200)[:5]
        nfs.append(f)
    nfs = np.array(nfs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    # adjust manually
#    j = 30
#    i,i2 = 1,3
#    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
#                       nfs[i,j+1:].copy())

    styles = ['-',':','--','.-','-.']
    labels = 'mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5'
    fig = figure(figsize=(4,3.6))
    ax = fig.add_subplot(111)
    for i in range(4):
        ax.loglog(ks,nfs[i],styles[i],label=labels[i])
    ax.set_xlabel(r'$k_{u2}\ (\mathrm{N/m})$', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
    ax.legend(loc='best')
#    xticks(fontsize=15)
#    yticks(fontsize=15)
    fig.tight_layout()
#    ax.set_ylim(0,300)
    savefig('example2_f_ku.png',dpi=800)
    savefig('example2_f_ku.svg',dpi=800)
    return nfs

def f_kv():
    ks = logspace(2,10,50)
    nfs = []
    for k in ks:
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam()
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].free()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].add_spring(0,k,0)
        f = beam.natural_frequency_real(0,500,200)[:5]
        nfs.append(f)
    nfs = np.array(nfs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    # adjust manually
#    j = 30
#    i,i2 = 1,3
#    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
#                       nfs[i,j+1:].copy())

    styles = ['-',':','--','.-','-.']
    labels = 'mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5'
    fig = figure(figsize=(4,3.6))
    ax = fig.add_subplot(111)
    for i in range(4):
        ax.loglog(ks,nfs[i],styles[i],label=labels[i])
    ax.set_xlabel(r'$k_{v2}\ (\mathrm{N/m})$', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
    ax.legend(loc='best')
#    xticks(fontsize=15)
#    yticks(fontsize=15)
    fig.tight_layout()
#    ax.set_ylim(0,300)
    savefig('example2_f_kv.png',dpi=800)
    savefig('example2_f_kv.svg',dpi=1000)
    return nfs


def f_kr():
    ks = logspace(2,10,50)
    nfs = []
    for k in ks:
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam()
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].free()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].add_spring(0,0,k)
        f = beam.natural_frequency_real(0,500,200)[:5]
        nfs.append(f)
    nfs = np.array(nfs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    # adjust manually
#    j = 30
#    i,i2 = 1,3
#    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
#                       nfs[i,j+1:].copy())

    styles = ['-',':','--','.-','-.']
    labels = 'mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5'
    fig = figure(figsize=(4,3.6))
    ax = fig.add_subplot(111)
    for i in range(4):
        ax.loglog(ks,nfs[i],styles[i],label=labels[i])
    ax.set_xlabel(r'$k_{\theta 2}\ (\mathrm{N/m})$', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
    ax.legend(loc='best')
#    xticks(fontsize=15)
#    yticks(fontsize=15)
    fig.tight_layout()
#    ax.set_ylim(0,300)
    savefig('example2_f_kr.png',dpi=800)
    savefig('example2_f_kr.svg',dpi=1000)
    return nfs




def forced_change_ku():
    Ea = 69e9 * (1+0.01j)
    Ec = 110e9 * (1+0.01j)
    sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                    [rhoc,rhoa], [0.85,0.85], 'section-1')
    sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                    [rhoc,rhoa], [0.85,0.85], 'section-2')
    sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                    [rhoc,rhoa], [0.85,0.85], 'section-3')
    beam = Beam()
    beam.append_segment(1.0, sect1)
    beam.append_segment(2.0, sect2)
    beam.append_segment(3.0, sect3)
    beam.nodes[0].add_spring(5e7,2e8)
    beam.nodes[1].add_mass(50)
    beam.nodes[2].add_mass(50)
    beam.nodes[3].add_spring(5e7,2e8)

    ks = [2e7,5e7,1e8]
    label_k = r'2\times 10^7', r'5\times 10^7', r'1\times 10^{8}'
    styles = ['-.','-','--']
    fig1, ax1 = subplots(figsize=(8,4))
    fig2, ax2 = subplots(figsize=(8,4))
    for i,k in enumerate(ks):
        beam.nodes[0].free().add_spring(k,2e8)
        beam.nodes[3].free().add_spring(k,2e8)
        freq = linspace(1,1000,1000)
        h1 = beam.transfer_function(1,'fx',2,'u',freq)
        h2 = beam.transfer_function(1,'fx',2,'v',freq)
        ax1.semilogy(freq,abs(h1),styles[i], label=r'$k_{{u1}}=k_{{u2}}={k}\mathrm{{N/m}}$'.format(k=label_k[i]))
        ax2.semilogy(freq,abs(h2),styles[i], label=r'$k_{{u1}}=k_{{u2}}={k}\mathrm{{N/m}}$'.format(k=label_k[i]))

    for ax in [ax1, ax2]:
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Displacement (m)')
        ax.legend(loc='best',fontsize=12)
        ax.set_xlim(0,1000)
    fig1.tight_layout()
    fig1.savefig('example2_forced_change_ku_u.png', dpi=800)
    fig1.savefig('example2_forced_change_ku_u.svg', dpi=1000)
    fig2.tight_layout()
    fig2.savefig('example2_forced_change_ku_v.png', dpi=800)
    fig2.savefig('example2_forced_change_ku_v.svg', dpi=1000)

def forced_change_kv():
    Ea = 69e9 * (1+0.01j)
    Ec = 110e9 * (1+0.01j)
    sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                    [rhoc,rhoa], [0.85,0.85], 'section-1')
    sect2 = Section([-0.1,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                    [rhoc,rhoa], [0.85,0.85], 'section-2')
    sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                    [rhoc,rhoa], [0.85,0.85], 'section-3')
    beam = Beam()
    beam.append_segment(1.0, sect1)
    beam.append_segment(2.0, sect2)
    beam.append_segment(3.0, sect3)
    beam.nodes[0].add_spring(5e7,2e8)
    beam.nodes[1].add_mass(50)
    beam.nodes[2].add_mass(50)
    beam.nodes[3].add_spring(5e7,2e8)

    ks = [5e7,1e8,2e8]
    label_k = r'5\times 10^7', r'1\times 10^8', r'2\times 10^{8}'
    styles = ['-.','-','--']
    fig1, ax1 = subplots(figsize=(8,4))
    fig2, ax2 = subplots(figsize=(8,4))
    for i,k in enumerate(ks):
        beam.nodes[0].free().add_spring(5e7,k)
        beam.nodes[3].free().add_spring(5e7,k)
        freq = linspace(1,1000,1000)
        h1 = beam.transfer_function(1,'fx',2,'u',freq)
        h2 = beam.transfer_function(1,'fx',2,'v',freq)
        ax1.semilogy(freq,abs(h1),styles[i], label=r'$k_{{v1}}=k_{{v2}}={k}\mathrm{{N/m}}$'.format(k=label_k[i]))
        ax2.semilogy(freq,abs(h2),styles[i], label=r'$k_{{v1}}=k_{{v2}}={k}\mathrm{{N/m}}$'.format(k=label_k[i]))

    for ax in [ax1, ax2]:
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Displacement (m)')
        ax.legend(loc='best',fontsize=12)
        ax.set_xlim(0,1000)
    fig1.tight_layout()
    fig1.savefig('example2_forced_change_kv_u.png', dpi=800)
    fig1.savefig('example2_forced_change_kv_u.svg', dpi=1000)
    fig2.tight_layout()
    fig2.savefig('example2_forced_change_kv_v.png', dpi=800)
    fig2.savefig('example2_forced_change_kv_v.svg', dpi=1000)



#if __name__ == '__main__':
#    from multiprocessing import Pool
#    p = Pool(processes = 3)
#    p.apply_async(c_ku)
#    p.apply_async(c_kv)
#    p.apply_async(c_kr)
#    p.apply_async(s_ku)
#    p.apply_async(s_kv)
#    p.apply_async(s_kr)
#    p.apply_async(f_ku)
#    p.apply_async(f_kv)
#    p.apply_async(f_kr)
#    p.apply_async(forced_change_ku)
#    p.apply_async(forced_change_kv)
#    p.close()
#    p.join()
#    show()



