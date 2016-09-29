import numpy as np
import matplotlib.pyplot as plt
import corr_funcs as cf

def test_epdecorrelation(v2true = 0.1,
                         delta = 0.3,
                         Nevent = 5000,
                         Npart = 100):
    '''
    Test event plane decorrelation
    v2true := true v2 which the samples will be generated with
    delta  := event plane decorrelation angle
    Nevent := Number of events generated
    Npart  := Numper of particles per region
    '''

    # Print running conditions
    print("--------------------------------------------")
    print(" v2true        : {}".format(v2true))
    print(" delta         : {}".format(delta))
    print(" Nevent        : {}".format(Nevent))
    print(" Npart         : {}".format(Npart))
    print("--------------------------------------------")

    # Generate random event planes in the north (N), mid (M), and south (S)
    Psi_M = np.random.uniform(-np.pi, np.pi, Nevent)
    Psi_N = Psi_M + delta
    Psi_S = Psi_M - delta

    # print('Psi_M: {}'.format(Psi_M))
    # print('Psi_N: {}'.format(Psi_N))
    # print('Psi_S: {}'.format(Psi_S))


    # Generate particle's for each event
    part_M = [cf.gen_part(v2true, Psi_M[i], 2, Npart) for i in range(Nevent)]
    part_N1 = [cf.gen_part(v2true, Psi_N[i], 2, Npart) for i in range(Nevent)]
    part_N2 = [cf.gen_part(v2true, Psi_N[i], 2, Npart) for i in range(Nevent)]
    part_S1 = [cf.gen_part(v2true, Psi_S[i], 2, Npart) for i in range(Nevent)]
    part_S2 = [cf.gen_part(v2true, Psi_S[i], 2, Npart) for i in range(Nevent)]

    # print('part_M:\n{}'.format(part_M))

    # Calculate the delta phi correlations for different combinations
    dphi_MN1 = np.concatenate([cf.calc_dphi(part_M[i], part_N1[i]) for i in range(Nevent)])
    dphi_MN2 = np.concatenate([cf.calc_dphi(part_M[i], part_N2[i]) for i in range(Nevent)])
    dphi_MS1 = np.concatenate([cf.calc_dphi(part_M[i], part_S1[i]) for i in range(Nevent)])
    dphi_MS2 = np.concatenate([cf.calc_dphi(part_M[i], part_S2[i]) for i in range(Nevent)])
    dphi_N1N2 = np.concatenate([cf.calc_dphi(part_N1[i], part_N2[i]) for i in range(Nevent)])
    dphi_S1S2 = np.concatenate([cf.calc_dphi(part_S1[i], part_S2[i]) for i in range(Nevent)])
    dphi_N1S1 = np.concatenate([cf.calc_dphi(part_N1[i], part_S1[i]) for i in range(Nevent)])

    # Calculate the c2 coefficients for each correlation
    c2_MN1 = cf.calc_v2(dphi_MN1)
    c2_MN2 = cf.calc_v2(dphi_MN2)
    c2_MS1 = cf.calc_v2(dphi_MS1)
    c2_MS2 = cf.calc_v2(dphi_MS2)
    c2_N1N2 = cf.calc_v2(dphi_N1N2)
    c2_S1S2 = cf.calc_v2(dphi_S1S2)
    c2_N1S1 = cf.calc_v2(dphi_N1S1)

    print('\n-- c2 values from 2 particle correlations--')
    print('c2_MN1: {}'.format(c2_MN1))
    print('c2_MN2: {}'.format(c2_MN2))
    print('c2_MS1: {}'.format(c2_MS1))
    print('c2_MS2: {}'.format(c2_MS2))
    print('c2_N1N2: {}'.format(c2_N1N2))
    print('c2_S1S2: {}'.format(c2_S1S2))
    print('c2_N1S1: {}'.format(c2_N1S1))


    # Calculate v2
    v2_MN1N2 = cf.calc_vn3sub(c2_MN1, c2_MN2, c2_N1N2)
    v2_MS1S2 = cf.calc_vn3sub(c2_MS1, c2_MS2, c2_S1S2)
    v2_MN1S1 = cf.calc_vn3sub(c2_MN1, c2_MS1, c2_N1S1)

    print('\n-- v2 values from 2 particle correlations--')
    print('v2_MN1N2: {}'.format(v2_MN1N2))
    print('v2_MS1S2: {}'.format(v2_MS1S2))
    print('v2_MN1S1: {}'.format(v2_MN1S1))




    # Make some plots
    fig, ax = plt.subplots(figsize=(6, 6))

    x = [0, 1, 2]
    y = np.vstack((v2_MN1N2, v2_MS1S2, v2_MN1S1))
    labels = ['MN1N2', 'MS1S2', 'MN1S1']
    # plt.xticks(x, labels, rotation='vertical')
    plt.xticks(x, labels, rotation=45)

    ax.set_ylim(0.05, 0.15)
    ax.set_xlim(-0.5, 2.5)

    ax.errorbar(x, y[:, 0], yerr=y[:, 1],
                ls='*', marker='o', ms=6, color='dodgerblue')

    ax.axhline(y=v2true, linestyle='--', color='crimson')


    ax.text(0.07, 0.95, r'$N_{event}=$' + '{}'.format(Nevent),
                   fontsize=10, transform=ax.transAxes)
    ax.text(0.07, 0.90, r'$N_{part}=$' + '{}'.format(Npart),
                   fontsize=10, transform=ax.transAxes)
    ax.text(0.07, 0.85, r'$\Delta=$' + '{:.2}'.format(delta),
                   fontsize=10, transform=ax.transAxes)
    ax.text(0.07, 0.80, r'$v_2^{true}=$' + '{:.2}'.format(v2true), color='crimson',
                   fontsize=10, transform=ax.transAxes)


    fig.savefig('pdfs/v2.pdf')
    plt.close(fig)



if __name__ == '__main__':

    # run with default values
    test_epdecorrelation()