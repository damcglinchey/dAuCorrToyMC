import numpy as np
import matplotlib.pyplot as plt
import corr_funcs as cf

def test_epdecorrelation(Nevent = 5000,
                         eta = np.array((-3.25, -2.0, 0.0, 2.0, 3.25)),
                         ldet = ['BBCS', 'FVTXS', 'CNT', 'FVTXN', 'BBCN'],
                         # Npart = np.array((10, 10, 10, 10, 10)),
                         Npart = np.array((74, 17, 5, 9, 10)),
                         v2true = np.array((0.1, 0.1, 0.1, 0.1, 0.1)),
                         sigd = 0.0001,
                         ):
    '''
    Test event plane decorrelation
    The code is set up to look in 5 separate regions defined by the
    eta values, labels, and Npart arrays input as arguments
    Nevent := Number of events generated
    eta    := Array of "detector" eta values
    Npart  := Array of the numper of particles per region
    v2true := Array of true v2 in each region
    sigd   := Sigma value for the Gaussian EP decorrelation 
    '''

    # Get the number of detectors
    ndet = len(eta)

    # Print running conditions
    print("--------------------------------------------")
    print(" Nevent        : {}".format(Nevent))
    print(" ndet          : {}".format(ndet))
    print(" ldet          : {}".format(ldet))
    print(" eta           : {}".format(eta))
    print(" Npart         : {}".format(Npart))
    print(" v2true        : {}".format(v2true))
    print(" sigd          : {}".format(sigd))
    print("--------------------------------------------")

    # Generate random event planes 
    # Assume linear decorrelation with random slope:
    # Psi(eta, sigd) = G(0, sigd)*eta + Psi0
    Psi0 = np.random.uniform(-np.pi, np.pi, Nevent)
    m = np.random.normal(0, sigd, len(Psi0))
    Psi = [(m*eta[i] + Psi0) for i in range(ndet)]

    # print('Psi0:\n{}'.format(Psi0))
    # print('m:\n{}'.format(m))
    # print('Psi:\n{}'.format(Psi))

    # Generate particle's for each event
    part = []
    for i in range(ndet):
        part.append([cf.gen_part(v2true[i], Psi[i][j], 2, Npart[i]) for j in range(Nevent)])

    # for i in range(ndet):
    #     print('Npart[{}] (Nevents: {}):\n{}'.format(i, 
    #                                                 len(part[i]), 
    #                                                 [len(part[i][j]) for j in range(Nevent)]))

    # Calculate the delta phi correlations for all combinations
    dphi = []
    for i in range(ndet):
        dphii = []
        for j in range(ndet):
            dphii.append(np.concatenate([cf.calc_dphi(part[i][k], part[j][k]) for k in range(Nevent)]))
        dphi.append(dphii)

    # print('dphi:\n{}'.format(dphi))



    # Calculate the c2 coefficients for each correlation
    c2 = []
    for i in range(ndet):
        c2tmp = []
        for j in range(ndet):
            c2tmp.append(cf.calc_v2(dphi[i][j]))
        c2.append(c2tmp)

    print('\n-- c2 values --')
    for i in range(ndet):
        for j in range(ndet):
            print('c2[{}][{}]({}-{})={}'.format(i, j, 
                                                ldet[i], ldet[j], 
                                                c2[i][j]))


    # Calculate v2
    v2_210 = cf.calc_vn3sub(c2[2][0], c2[2][1], c2[0][1])
    v2_234 = cf.calc_vn3sub(c2[2][3], c2[2][4], c2[3][4])
    v2_213 = cf.calc_vn3sub(c2[2][1], c2[2][3], c2[1][3])
    v2_204 = cf.calc_vn3sub(c2[2][0], c2[2][4], c2[0][4])

    print('\n-- v2 values --')
    print('v2({}-{}-{}) = {}'.format(ldet[2], ldet[1], ldet[0], v2_210))
    print('v2({}-{}-{}) = {}'.format(ldet[2], ldet[3], ldet[4], v2_234))
    print('v2({}-{}-{}) = {}'.format(ldet[2], ldet[1], ldet[3], v2_213))
    print('v2({}-{}-{}) = {}'.format(ldet[2], ldet[0], ldet[4], v2_204))




    # Make some plots
    fig, ax = plt.subplots(figsize=(6, 6))

    x = [0, 1, 2, 3]
    y = np.vstack((v2_210, v2_234, v2_213, v2_204))
    labels = ['{}-{}-{}'.format(ldet[2], ldet[1], ldet[0]),
              '{}-{}-{}'.format(ldet[2], ldet[3], ldet[4]),
              '{}-{}-{}'.format(ldet[2], ldet[1], ldet[3]),
              '{}-{}-{}'.format(ldet[2], ldet[0], ldet[4])]
    # plt.xticks(x, labels, rotation='vertical')
    plt.xticks(x, labels, rotation=60, size=5)

    ax.set_ylim(0.05, 0.15)
    ax.set_xlim(-0.5, 3.5)

    ax.errorbar(x, y[:, 0], yerr=y[:, 1],
                ls='*', marker='o', ms=6, color='dodgerblue')

    ax.axhline(y=v2true[2], linestyle='--', color='crimson')


    ax.text(0.07, 0.95, r'$N_{event}=$' + '{}'.format(Nevent),
                   fontsize=10, transform=ax.transAxes)
    ax.text(0.07, 0.90, r'$\sigma_{\Delta}=$' + '{:.2}'.format(sigd),
                   fontsize=10, transform=ax.transAxes)
    ax.text(0.07, 0.85, r'$v_2^{true}=$' + '{:.2}'.format(v2true[2]), color='crimson',
                   fontsize=10, transform=ax.transAxes)


    fig.savefig('pdfs/v2.pdf')
    plt.close(fig)



if __name__ == '__main__':

    # run with default values
    test_epdecorrelation()