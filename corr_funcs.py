import numpy as np
import matplotlib.pyplot as plt


def gen_part(vn, psi, order=2, N=100):
    '''
    Generate random set of particles
    '''
    rand = np.zeros(N)
    for i in range(N):
        acc = False
        while acc is False:
            rphi = np.random.uniform(-np.pi, np.pi, 1)
            vr = 1 + 2*vn*np.cos(order*(rphi-psi))
            rv = np.random.uniform(0, 1+2*vn, 1)
            if rv < vr:
                rand[i] = rphi
                acc = True
    return rand


def gen_psi(Psi0, eta, sigd):
    '''
    Generate the event plane angles
    Assume linear decorrelation with random slope:
    Psi(eta, sigd) = G(0, sigd)*eta + Psi0
    '''
    # Generate random slope
    m = np.random.normal(0, sigd, len(Psi0))
    # Calculate psi
    Psi = m*eta + Psi0
    return Psi


def calc_cn(d, order):
    '''
    Calculate the correlation coefficients numerically
    From the data.
    d[:, 0] := array of delta phi values
    d[:, 1] := array of correlation values
    d[:, 2] := statistical uncertainties
    '''
    val = np.sum(d[:, 1] * np.cos(order*d[:, 0])) / np.sum(d[:, 1])
    err = np.sqrt(np.sum((d[:, 2] * np.cos(order*d[:, 0]))**2)) / np.sum(d[:, 1])
    return np.array((val, err))

def calc_v2(dphi):
    '''
    Calculate the v2 (or c2) value from delta phi values along
    with it's uncertainty.
    dphi := array of delta phi values
    '''
    v2 = np.mean(np.cos(2*(dphi)))
    v2e = np.sqrt(np.sum((v2 - np.cos(2*dphi))**2)) / len(dphi)
    return np.array((v2, v2e))

def calc_dphi(phi1, phi2):
    '''
    Calculate delta phi values for all combinations of particles
    phi1 := array of phi_1 values
    phi2 := array of phi_2 values
    '''
    dphi = np.zeros(len(phi1)*len(phi2))
    i = 0
    for p1 in phi1:
        for p2 in phi2:
            dp = p1 - p2
            if dp < -np.pi:
                dp = dp + np.pi
            if dp > np.pi:
                dp = dp - np.pi
            dphi[i] = dp
            i = i + 1
    return dphi

def calc_vn3sub(cAB, cAC, cBC):
    ''' 
    Calculate the vn using 3 sub-event method from 2 particle correlations
    '''
    val = np.sqrt((cAB[0]*cAC[0])/cBC[0])

    # calculate the error
    dcAB = 0.5 * (cAC[0] / cBC[0]) / val
    dcAC = 0.5 * (cAB[0] / cBC[0]) / val
    dcBC = -0.5 * (cAB[0] * cAC[0]) / (cBC[0]**2) / val

    err = np.sqrt((dcAB*cAB[1])**2 + (dcAC*cAC[1])**2 + (dcBC*cBC[1])**2)

    return np.array((val, err))





if __name__ == '__main__':

    # test generation
    vn = 0.1
    psi = np.pi / 2.
    phi = np.arange(-np.pi, np.pi, 2*np.pi/100.)
    Ndphi = 100*(1 + 2*vn*np.cos(2*(phi - psi)))

    # random phi values
    rphi = gen_part(vn, psi, 2, 10000)

    # calculate v2 from random phi values
    v2obs = calc_v2(rphi - psi)


    fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot(phi, Ndphi, ls='--', color='dodgerblue')
    ax.hist(rphi, 100, histtype='stepfilled', facecolor='gray', alpha=0.6)

    ax.text(0.07, 0.95, r'$v_2^{true}=$' + '{:.2}'.format(vn),
                   fontsize=6, transform=ax.transAxes)
    obs = r'$v_2^{obs}=$' + '{:.2}'.format(v2obs[0]) + r'$\pm$' + '{:.2}'.format(v2obs[1])
    ax.text(0.07, 0.90, obs, fontsize=6, transform=ax.transAxes)

    fig.savefig('test_gen.pdf')
    plt.close(fig)

