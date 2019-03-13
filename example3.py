#!/usr/bin/env python3


from pylab import *
from beam import *
import matplotlib.ticker as ticker

def change_hc2_cc():
    fs = []
    hs = linspace(0,0.3,61)
    for h in hs:
        Es = 210e9
        Ea = 69e9
        rhos = 7800
        rhoa = 2700
        nus = 0.3
        nua = 0.33
        Ec = 110e9
        nuc = 0.32
        rhoc = 8500
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-h,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam(P=0)
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].fixed()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].fixed()
        f = beam.natural_frequency_real(0,800,800)[:6]
        fs.append(f)
    nfs = np.array(fs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
    #                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
#    j = 28 #28
#    i, i2 = 3, 5
#    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
#                       nfs[i,j+1:].copy())
    for i in range(5):
        fig,ax=subplots(figsize=(3,3))
        ax.plot(hs, nfs[i])
        ax.set_xlabel(r'$h_{c2}\ \mathrm{(m)}$', fontsize=12)
        ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
        fig.tight_layout()
        savefig('example3_cc_{n}.png'.format(n=i+1), dpi=800)
        savefig('example3_cc_{n}.svg'.format(n=i+1), dpi=1000)

def change_hc2_ss():
    fs = []
    hs = linspace(0,0.3,61)
    for h in hs:
        Es = 210e9
        Ea = 69e9
        rhos = 7800
        rhoa = 2700
        nus = 0.3
        nua = 0.33
        Ec = 110e9
        nuc = 0.32
        rhoc = 8500
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-h,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam(P=0)
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].pinned()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].pinned()
        f = beam.natural_frequency_real(0,800,800)[:6]
        fs.append(f)
    nfs = np.array(fs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    for i in range(5):
        fig,ax=subplots(figsize=(3,3))
        ax.plot(hs, nfs[i])
        ax.set_xlabel(r'$h_{c2}\ \mathrm{(m)}$', fontsize=12)
        ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
        fig.tight_layout()
        savefig('example3_ss_{n}.png'.format(n=i+1), dpi=800)
        savefig('example3_ss_{n}.svg'.format(n=i+1), dpi=1000)

def change_hc2_ff():
    fs = []
    hs = linspace(0,0.3,61)
    for h in hs:
        Es = 210e9
        Ea = 69e9
        rhos = 7800
        rhoa = 2700
        nus = 0.3
        nua = 0.33
        Ec = 110e9
        nuc = 0.32
        rhoc = 8500
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-h,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam(P=0)
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].free()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].free()
        f = beam.natural_frequency_real(0,800,800)[:6]
        fs.append(f)
    nfs = np.array(fs).T
    for j in range(1,nfs.shape[1]-1):
        for i in range(len(nfs)):
            for i2 in range(len(nfs)):
                if (abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i,j+1]) >
                        abs(2*nfs[i,j] - nfs[i,j-1] -nfs[i2,j+1])):
#                    print(j,i,i2)
                    nfs[i,j+1:],nfs[i2,j+1:] = (nfs[i2,j+1:].copy(),
                       nfs[i,j+1:].copy())
    for i in range(5):
        fig,ax=subplots(figsize=(3,3))
        ax.plot(hs, nfs[i])
        ax.set_xlabel(r'$h_{c2}\ \mathrm{(m)}$', fontsize=12)
        ax.set_ylabel(r'$\mathrm{Natural\ frequency\ (Hz)}$', fontsize=12)
        fig.tight_layout()
        savefig('example3_ff_{n}.png'.format(n=i+1), dpi=800)
        savefig('example3_ff_{n}.svg'.format(n=i+1), dpi=1000)


def buckling_find_p_range(beam, p_min, p_max, Ps, ffs, n=50):
#    print(p_min, p_max)
    EPS = 0
    if len(Ps) > n:
        return Ps, ffs
    p = linspace(p_min, p_max, n//2+1, endpoint=False)
    for i in range(1, len(p)):
        beam.reset_P([-p[i]]*3)
        f = beam.natural_frequency_real(0,ffs[-1]+EPS,50)
        if len(f) == 0:
            p_max = p[i]
            p_min = p[i-1]
            return buckling_find_p_range(beam, p_min, p_max, Ps, ffs, n)
        f = f[0]
        Ps.append(p[i])
        ffs.append(f)
    p_min = p[-1]
    return buckling_find_p_range(beam, p_min, p_max, Ps, ffs, n)

def buckling_cc():
    Es = 210e9
    Ea = 69e9
    rhos = 7800
    rhoa = 2700
    nus = 0.3
    nua = 0.33
    Ec = 110e9
    nuc = 0.32
    rhoc = 8500
    hs = [0.0,0.1,0.2,0.3]
    fig, ax = subplots(figsize=(3,3))
    styles = '-','-.','--',':'
    for i,h in enumerate(hs):
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-h,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam(P=0)
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].fixed()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].fixed()
        P_min = 0
        P_max = [2e6,6e6,8e6,1e7][i]
        Ps = [P_min]
        beam.reset_P([-P_min,-P_min,-P_min])
        ffs = [beam.natural_frequency_real(0,100,20)[0]]
        Ps, ffs = buckling_find_p_range(beam, P_min, P_max, Ps, ffs,50)
        ax.plot(Ps, ffs, styles[i],label=r'$h_{{c2}}={h}\mathrm{{m}}$'.format(h=h))
    ax.set_xlabel(r'$\mathrm{Axial\ load\ (N)}$')
    ax.set_ylabel(r'$\mathrm{Fundamental\ frequency\ (Hz)}$')
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.legend(loc='best')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
    fig.tight_layout()
    savefig('example3_buckling_cc.png',dpi=800)
    savefig('example3_buckling_cc.svg',dpi=1000)


def buckling_ss():
    Es = 210e9
    Ea = 69e9
    rhos = 7800
    rhoa = 2700
    nus = 0.3
    nua = 0.33
    Ec = 110e9
    nuc = 0.32
    rhoc = 8500
    hs = [0.0,0.1,0.2,0.3]
    fig, ax = subplots(figsize=(3,3))
    styles = '-','-.','--',':'
    for i,h in enumerate(hs):
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-h,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam(P=0)
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].pinned()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].pinned()
        P_min = 0
        P_max = [2e6,6e6,8e6,1e7][i]
        Ps = [P_min]
        beam.reset_P([-P_min,-P_min,-P_min])
        ffs = [beam.natural_frequency_real(0,100,20)[0]]
        Ps, ffs = buckling_find_p_range(beam, P_min, P_max, Ps, ffs,50)
        ax.plot(Ps, ffs, styles[i],label=r'$h_{{c2}}={h}\mathrm{{m}}$'.format(h=h))
    ax.set_xlabel(r'$\mathrm{Axial\ load\ (N)}$')
    ax.set_ylabel(r'$\mathrm{Fundamental\ frequency\ (Hz)}$')
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.legend(loc='best')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
    fig.tight_layout()
    savefig('example3_buckling_ss.png',dpi=800)
    savefig('example3_buckling_ss.svg',dpi=1000)

def buckling_ff():
    Es = 210e9
    Ea = 69e9
    rhos = 7800
    rhoa = 2700
    nus = 0.3
    nua = 0.33
    Ec = 110e9
    nuc = 0.32
    rhoc = 8500
    hs = [0.0,0.1,0.2,0.3]
    fig, ax = subplots(figsize=(3,3))
    styles = '-','-.','--',':'
    for i,h in enumerate(hs):
        sect1 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-1')
        sect2 = Section([-h,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-2')
        sect3 = Section([-0.05,0,0.05], [Ec,Ea], [nuc, nua], [0.1,0.1],
                        [rhoc,rhoa], [0.85,0.85], 'section-3')
        beam = Beam(P=0)
        beam.append_segment(1.0, sect1)
        beam.append_segment(2.0, sect2)
        beam.append_segment(3.0, sect3)
        beam.nodes[0].free()
        beam.nodes[1].add_mass(50)
        beam.nodes[2].add_mass(50)
        beam.nodes[3].free()
        P_min = 0
        P_max = [2e6,6e6,8e6,1e7][i]
        Ps = [P_min]
        beam.reset_P([-P_min,-P_min,-P_min])
        ffs = [beam.natural_frequency_real(0,100,20)[0]]
        Ps, ffs = buckling_find_p_range(beam, P_min, P_max, Ps, ffs,50)
        ax.plot(Ps, ffs, styles[i],label=r'$h_{{c2}}={h}\mathrm{{m}}$'.format(h=h))
    ax.set_xlabel(r'$\mathrm{Axial\ load\ (N)}$')
    ax.set_ylabel(r'$\mathrm{Fundamental\ frequency\ (Hz)}$')
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.legend(loc='best')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
    fig.tight_layout()
    savefig('example3_buckling_ff.png',dpi=800)
    savefig('example3_buckling_ff.svg',dpi=1000)







#if __name__ == '__main__':
#    from multiprocessing import Pool
#    p = Pool(processes = 3)
###    p.apply_async(change_hc2_cc)
###    p.apply_async(change_hc2_ss)
###    p.apply_async(change_hc2_ff)
#    p.apply_async(buckling_cc)
#    p.apply_async(buckling_ss)
#    p.apply_async(buckling_ff)
#
#    p.close()
#    p.join()
#    show()
