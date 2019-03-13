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
#
#f = beam.natural_frequency_real(0,2000,5000)

def natural_frequency():
    f = beam.natural_frequency_real(0,2000)
    print('固有频率为：')
    print(f)


def mode_shape():
    m = beam.modal(0,2000)
    for i in range(10):
        fig = plotbeam(m[i],figsize=(6,3),scalar=0.3)
        ax = fig.axes[0]
        ax.set_ylim(-0.4,0.4)
        fig.tight_layout()
        savefig('example1_mode_shape_{n}.png'.format(n=i+1), dpi=800)
        savefig('example1_mode_shape_{n}.svg'.format(n=i+1), dpi=1000)



def read_abaqus_result(filename='../abaqus/example 1/abaqus.rpt'):
    freq = []
    u1 = []
    u2 = []
    with open(filename) as f:
        for line in f:
            words = line.split()
            if len(words) != 3:
                continue
            try:
                words = [eval(i) for i in words][:3]
            except Exception:
                continue
            freq.append(words[0])
            u1.append(words[1])
            u2.append(words[2])
    res = np.array([freq,u1,u2])
    return res.T


def harmonic_response():
    Es = 210e9 * (1+0.01j)
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

    res_abaqus = read_abaqus_result()
    freq = linspace(1,1000,1000)
    h1 = beam.transfer_function(1,'fx',2,'u',freq)
    h2 = beam.transfer_function(1,'fx',2,'v',freq)

    fig, ax = subplots(figsize=(8,4))
    ax.semilogy(freq,res_abaqus[:,1],'-',lw=5, label='FEM')
    ax.semilogy(freq,abs(h1),'k-', label='Presented method')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Displacement (m)')
    ax.legend(loc='best',fontsize=12)
    ax.set_xlim(0,1000)
    fig.tight_layout()
    savefig('example1_harmonic_response_u.png', dpi=800)
    savefig('example1_harmonic_response_u.svg', dpi=1000)

    fig, ax = subplots(figsize=(8,4))
    ax.semilogy(freq,res_abaqus[:,2],'-',lw=5, label='FEM')
    ax.semilogy(freq,abs(h2),'k-', label='Presented method')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Displacement (m)')
    ax.legend(loc='best',fontsize=12)
    ax.set_xlim(0,1000)
    fig.tight_layout()
    savefig('example1_harmonic_response_v.png', dpi=800)
    savefig('example1_harmonic_response_v.svg', dpi=1000)

def harmonic_response_detail():
    Es = 210e9 * (1+0.01j)
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

    freqs = [17,60,89,197,273,402,654,683,720,980]
    for freq in freqs:
        beam.set_frequency(freq)
        res = beam.harmonic_response(1,1)
        fig = plotbeam(res, figsize=(6,2))
        ax = fig.axes[0]
        ax.axis('off')
        ax.set_ylim(-.4,.4)
        fig.tight_layout()
        savefig('example1_harmonic_response_{f}Hz.png'.format(f=freq), dpi=800)
        savefig('example1_harmonic_response_{f}Hz.svg'.format(f=freq), dpi=1000)

#natural_frequency()
#mode_shape()
harmonic_response()
#harmonic_response_detail()

show()









