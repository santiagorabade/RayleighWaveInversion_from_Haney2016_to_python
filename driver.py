def raylee_lysmer(Nn, vsv, vpv, rhov, f, hv, modn, Nnf, vpfv, rhofv, hfv):
# the number of nodes in the fluid, based on the number of elements

    if Nnf > 0:
        Nnfo = Nnf + 1
    else:
        Nnfo = 0

    # make fluid portion of model

    # make kappaf, the modulus
    kappafv = rhofv * vpfv * vpfv

    # make angular frequency
    omga = 2 * np.pi * f

    # initialize some local matrices
    L1 = np.zeros((2, 2))
    L3 = np.zeros((2, 2))
    M1 = np.zeros((2, 2))
    # initialize the global matrix
    Ka1 = np.zeros((Nnfo + (2 * Nn), Nnfo + (2 * Nn)))
    Ka3 = np.zeros((Nnfo + (2 * Nn), Nnfo + (2 * Nn)))
    M = np.zeros((Nnfo + (2 * Nn), Nnfo + (2 * Nn)))
    # for all elements
    for ii in range(1, Nnf+1):
        # grab grid interval of current element
        h = hfv[ii-1]
        
        # grab material properties of current element
        rhof = rhofv[ii-1]
        kappaf = kappafv[ii-1]

        # make elemental mass matrix
        M1 = np.zeros((2, 2))
        M1[0, 0] = h / (2 * kappaf)
        M1[1, 1] = h / (2 * kappaf)

        # make elemental stiffness matrices
        L1 = np.zeros((2, 2))
        L3 = np.zeros((2, 2))

        # some alternate variables from Lysmer
        alph = 1 / (6 * rhof)
        bet = 1 / (6 * rhof)


        # the 4 entries of the 2x2 elemental stiffness matrices of Lysmer
        L1[0, 0] = 2 * alph * h
        L3[0, 0] = (6 * bet / h)
        L1[0, 1] = alph * h
        L3[0, 1] = -(6 * bet / h)
        L1[1, 0] = L1[0, 1]
        L3[1, 0] = L3[0, 1]
        L1[1, 1] = L1[0, 0]
        L3[1, 1] = L3[0, 0]

        # assemble mass and stiffness matrices from elemental matrices
        M[ii-1:ii+1, ii-1:ii+1] += M1
        Ka1[ii-1:ii+1, ii-1:ii+1] += L1
        Ka3[ii-1:ii+1, ii-1:ii+1] += L3

    # make solid portion of model
    M[0, 0] *= 2
    Ka1[0, 0] *= 2
    Ka3[0, 0] *= 2
    # make mu and lambda
    muv = rhov * vsv * vsv
    lamdav = rhov * vpv * vpv - 2 * muv

    # initialize some matrices
    Ka2 = np.zeros((Nnfo + (2 * Nn), Nnfo + (2 * Nn)))

    L1 = np.zeros((4, 4))
    L2 = np.zeros((4, 4))
    L3 = np.zeros((4, 4))
    M1 = np.zeros((4, 4))

    # for all elements
    for ii in range(Nn):
        # grab grid interval of current element
        h = hv[ii]

        # grab material properties of current element
        mu = muv[ii]
        lamda = lamdav[ii]
        # make elemental mass matrix
        M1 = np.zeros((4, 4))
        M1[0, 0] = h * rhov[ii] / 2
        M1[1, 1] = h * rhov[ii] / 2
        M1[2, 2] = h * rhov[ii] / 2
        M1[3, 3] = h * rhov[ii] / 2

        # make elemental stiffness matrices
        L1 = np.zeros((4, 4))
        L2 = np.zeros((4, 4))
        L3 = np.zeros((4, 4))

        # some alternate variables from Lysmer
        alph = ((2 * mu) + lamda) / 6
        bet = mu / 6
        theta = (mu + lamda) / 4
        psi = (mu - lamda) / 4

        # the 16 entries of the 4x4 elemental stiffness matrices of Lysmer
        L1[0, 0] = 2 * alph * h
        L3[0, 0] = (6 * bet / h)

        L2[0, 1] = 2 * psi

        L1[0, 2] = alph * h
        L3[0, 2] = -(6 * bet / h)

        L2[0, 3] = 2 * theta

        L2[1, 0] = L2[0, 1]

        L1[1, 1] = 2 * bet * h
        L3[1, 1] = (6 * alph / h)

        L2[1, 2] = -2 * theta

        L1[1, 3] = bet * h
        L3[1, 3] = -(6 * alph / h)

        L1[2, 0] = L1[0, 2]
        L3[2, 0] = L3[0, 2]

        L2[2, 1] = L2[1, 2]

        L1[2, 2] = L1[0, 0]
        L3[2, 2] = L3[0, 0]

        L2[2, 3] = -2 * psi

        L2[3, 0] = L2[0, 3]

        L1[3, 1] = L1[1, 3]
        L3[3, 1] = L3[1, 3]

        L2[3, 2] = L2[2, 3]

        L1[3, 3] = L1[1, 1]
        L3[3, 3] = L3[1, 1]
        if ii==Nn-1:
            M[Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+1,Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+1] += M1[0:2,0:2]
            Ka1[Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+1,Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+1] += L1[0:2,0:2]
            Ka2[Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+1,Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+1] += L2[0:2,0:2]
            Ka3[Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+1,Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+1] += L3[0:2,0:2]
        else:
            M[Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+2,Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+2] += M1
            Ka1[Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+2,Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+2] += L1
            Ka2[Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+2,Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+2] += L2
            Ka3[Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+2,Nnfo+(2*(ii)):Nnfo+(2*(ii+1))+2] += L3

    # Construct the coupling matrix
    if Nnf > 0:
        Cm = np.zeros((Nnfo + (2 * Nn), Nnfo + (2 * Nn)))
        Cm[Nnfo, Nnfo + 2] = 1
        Cm[Nnfo + 2, Nnfo] = 1
    else:
        Cm = np.zeros((Nnfo + (2 * Nn), Nnfo + (2 * Nn)))

    # Find the Rayleigh/Scholte wave speed
    if Nnf > 0:
        msval, msloc = np.min(vsv), np.argmin(vsv)
        mfval, mfloc = np.min(vpfv), np.argmin(vpfv)
        vsmay = msval
        vpmay = vpv[msloc]
        vpfmay = mfval
        rhofmay = rhofv[mfloc]
        rhomay = rhov[msloc]
        rspd = stoneley_vel(vpmay, vsmay, vpfmay, rhofmay, rhomay)
    else:
        msval, msloc = np.min(vsv), np.argmin(vsv)
        vsmay = msval
        vpmay = vpv[msloc]
        t1 = 1 / (vsmay ** 6)
        t2 = -8 / (vsmay ** 4)
        t3 = (24 / (vsmay ** 2)) - (16 / (vpmay ** 2))
        t4 = -16 * (1 - ((vsmay / vpmay) ** 2))
        rspd = np.sqrt(np.min(np.roots([t1, t2, t3, t4])))
    
    # Find the eigenvalue closest to the upper-bound eigenvalue
    mn = int(modn)
    block1 = np.hstack([np.zeros((Nnfo+(2*Nn), Nnfo+(2*Nn))), np.eye(Nnfo+(2*Nn))])
    block2 = np.hstack([(omga*omga*M) - Ka3 - (omga*Cm), Ka2])
    block1D = np.diagonal(block1)
    block3 = np.hstack([np.eye(Nnfo+(2*Nn)), np.zeros((Nnfo+(2*Nn), Nnfo+(2*Nn)))])
    block4 = np.hstack([np.zeros((Nnfo+(2*Nn), Nnfo+(2*Nn))), Ka1])    
    block3D = np.diagonal(block3)
    sigmax = omga/rspd
    blockT = np.vstack([block1, block2])
    blockB = np.vstack([block3, block4])
    dp, xp = eigs(np.vstack([block1, block2]),k=mn,
                  M=np.vstack([block3, block4]),
                  sigma=sigmax)
    if mn>1:
        x = np.real(xp[:,mn-1])
        d = np.array([np.real(dp[mn-1])])
    else:
        x = np.real(xp)
        d = np.real(dp)
    # Normalize the eigenfunction
    fctr = (1 / (np.transpose(x[0:(Nnfo+(2*Nn))]) @ M @ x[0:(Nnfo+(2*Nn))])) \
    - ((np.transpose(x[0:(Nnfo+(2*Nn))]) @ Cm @ x[0:(Nnfo+(2*Nn))]) / (2 * omga))

    evp = x[0:(Nnfo+(2*Nn))] * np.sqrt(fctr) * np.sign(x[Nnfo])
    ev = evp[Nnfo:(Nnfo+(2*Nn))]
    evAux = ev
    # Calculate the wavenumber
    kk = d

    # Calculate the phase velocity
    vpk = omga / kk
    # Calculate the group velocity
    vgk = (np.transpose(x[0:Nnfo+(2*Nn)]).dot((2*d[0]*Ka1) - Ka2).dot(x[0:Nnfo+(2*Nn)])) / \
          ((2*omga*(np.transpose(x[0:Nnfo+(2*Nn)]).dot(M).dot(x[0:Nnfo+(2*Nn)]))) - np.transpose(x[0:Nnfo+(2*Nn)]).dot(Cm).dot(x[0:Nnfo+(2*Nn)]))
    # Test if it's a guided mode
    a = np.abs(ev[1::2])
    hs = np.concatenate(([0], np.cumsum(hv)))
    srti = np.argsort(np.abs(np.concatenate(([0], np.cumsum(hv))) - (np.sum(hv) / 2)))
    srtv = np.sort(np.abs(np.concatenate(([0], np.cumsum(hv))) - (np.sum(hv) / 2)))
    s3 = np.min(srti[0:1])

    dum2 = np.sum(a[s3:Nn].reshape(-1)*hv[s3:Nn])
    dum3 = np.sum(a[0:s3].reshape(-1)*hv[0:s3])


    srtvv = np.sort(np.abs(hs - (modn * 0.5 * (vpk / f))))
    srtii = np.argsort(np.abs(np.concatenate(([0], np.cumsum(hv))) - (modn * 0.5 * (vpk / f))))
    srtvv = np.sort(np.abs(np.concatenate(([0], np.cumsum(hv))) - (modn * 0.5 * (vpk / f))))
    if dum3 / dum2 < 3:
        vpk = np.nan
        vgk = np.nan
        kk = np.nan
        ev = np.full((2 * Nn, 1), np.nan)
    elif (2 * modn * 0.5 * (vpk / f)) / np.sum(hv) > 1.6:
        print((2 * modn * 0.5 * (vpk / f)) / np.sum(hv))
        raise ValueError("Insufficiently deep grid: Base of model less than twice the sensitivity depth from the top. Extend grid in depth and re-run.")
    elif np.sum((vpk / f) / hv[0:srtii[1]]) < 3:
        print(np.sum((vpk / f) / hv[0:srtii[1]]))
        raise ValueError("Insufficiently dense grid: Number of elements per wavelength less than 3 above the depth of sensitivity. Densify grid and re-run.")
    elif Nnf > 0:
        if (vpk / f) / np.max(hfv) < 5:
            raise ValueError("Insufficiently dense grid: Number of elements per wavelength less than 5 in fluid. Densify grid and re-run.")
    else:
        # The mode is acceptable, a guided mode
        pass
    
    return kk, vpk, vgk, ev

def raylee_sensitivity(Nn, vsv, vpv, rhov, fks, hv, modnv, vflg, Nnf, vpfv, rhofv, hfv, pratioflag):
    countr = 0

    #####################################################################
    # augment the frequency vector if group kernels are needed, 
    # this triples the size of the frequency vector
    #####################################################################
    fks_orig = fks
    #modnv_orig = modnv
    modnv_orig = modn
    if vflg[0] == 1:
        fks = [fks_orig[0]*.999, fks_orig[0], fks_orig[0]*1.001]
        modnv = [modnv_orig[0], modnv_orig[0], modnv_orig[0]]
    else:
        fks = fks_orig[0]
        modnv = modnv_orig[0]
    for ii in range(1, len(fks_orig)):
        if vflg[ii] == 1:
            fks.extend([fks_orig[ii]*.999, fks_orig[ii], fks_orig[ii]*1.001])
            modnv.extend([modnv_orig[ii], modnv_orig[ii], modnv_orig[ii]])
        else:
            fks=np.append(fks,fks_orig[ii])
            modnv=np.append(modnv,modnv_orig[ii])
    #####################################################################
    # calculate phase velocities and eigenfunctions
    #####################################################################
    vp = np.zeros(len(fks))
    U = np.zeros(len(fks))
    xm = np.zeros((2 * Nn, len(fks)))
    dum2p = np.zeros(len(fks))
    dum3p = np.zeros(len(fks))
    for f in fks:
        result = raylee_lysmer(Nn, vsv, vpv, rhov, f, hv, modnv[countr], Nnf, vpfv, rhofv, hfv)
        kk, vpk, vgk, x = result
        vp[countr] = vpk
        U[countr] = vgk
        for i in range(len(x)):
            xm[i][countr] = x[i]
        countr += 1    
    # Construct phase sensitivity kernels
    # Initialize matrices with zeros
    snsmf = np.zeros((Nn, len(fks)))
    snsmfh = np.zeros((Nn, len(fks)))
    snsmflam = np.zeros((Nn, len(fks)))
    snsmfrho = np.zeros((Nn, len(fks)))

    # For all depths (except the bottom element) and frequencies
    for ii in range(Nn-2):

        h = hv[ii]

        countr = 0
        for f in fks:

            # Density sensitivity
            snsmfrho[ii, countr] = -(vp[countr] / (2 * U[countr])) * \
                                   (xm[2*ii, countr] * (h/2) * xm[2*ii, countr] +
                                    xm[2*ii+1, countr] * (h/2) * xm[2*ii+1, countr] +
                                    xm[2*ii+2, countr] * (h/2) * xm[2*ii+2, countr] +
                                    xm[2*ii+3, countr] * (h/2) * xm[2*ii+3, countr]) * \
                                   rhov[ii]

            # Lambda sensitivity
            # k^2 term
            snsmflam[ii, countr] = (1 / (2 * vp[countr] * U[countr])) * \
                                   (xm[2*ii, countr] * (h/3) * xm[2*ii, countr] + \
                                    xm[2*ii+2, countr] * (h/3) * xm[2*ii+2, countr] + \
                                    xm[2*ii+2, countr] * (h/6) * xm[2*ii, countr] + \
                                    xm[2*ii, countr] * (h/6) * xm[2*ii+2, countr]) * \
                                    (rhov[ii] * ((vpv[ii]**2) - (2 * (vsv[ii]**2))))
            # k^1 term
            snsmflam[ii, countr] = snsmflam[ii, countr] + \
                                   (1 / (2 * (2 * np.pi * f) * U[countr])) * \
                                   (xm[2*ii, countr] * (1/2) * xm[2*ii+1, countr] + \
                                    xm[2*ii, countr] * (-1/2) * xm[2*ii+3, countr] + \
                                    xm[2*ii+1, countr] * (1/2) * xm[2*ii, countr] + \
                                    xm[2*ii+1, countr] * (1/2) * xm[2*ii+2, countr] + \
                                    xm[2*ii+2, countr] * (1/2) * xm[2*ii+1, countr] + \
                                    xm[2*ii+2, countr] * (-1/2) * xm[2*ii+3, countr] + \
                                    xm[2*ii+3, countr] * (-1/2) * xm[2*ii, countr] + \
                                    xm[2*ii+3, countr] * (-1/2) * xm[2*ii+2, countr]) * \
                                    (rhov[ii] * ((vpv[ii]**2) - (2 * (vsv[ii]**2))))
            
            # k^0 term
            snsmflam[ii, countr] = snsmflam[ii, countr] + \
                                   (vp[countr] / (2 * ((2 * np.pi * f)**2) * U[countr])) * \
                                   (xm[2*ii+1, countr] * (1/h) * xm[2*ii+1, countr] + \
                                    xm[2*ii+3, countr] * (1/h) * xm[2*ii+3, countr] + \
                                    xm[2*ii+1, countr] * (-1/h) * xm[2*ii+3, countr] + \
                                    xm[2*ii+3, countr] * (-1/h) * xm[2*ii+1, countr]) * \
                                    (rhov[ii] * ((vpv[ii]**2) - (2 * (vsv[ii]**2))))


            # mu sensitivity
            # k^2 term
            snsmf[ii, countr] = (1 / (2 * vp[countr] * U[countr])) * \
                                (xm[2*ii, countr] * (2*h/3) * xm[2*ii, countr] + \
                                 xm[2*ii+2, countr] * (2*h/3) * xm[2*ii+2, countr] + \
                                 xm[2*ii, countr] * (h/3) * xm[2*ii+2, countr] + \
                                 xm[2*ii+2, countr] * (h/3) * xm[2*ii, countr] + \
                                 xm[2*ii+1, countr] * (h/3) * xm[2*ii+1, countr] + \
                                 xm[2*ii+3, countr] * (h/3) * xm[2*ii+3, countr] + \
                                 xm[2*ii+1, countr] * (h/6) * xm[2*ii+3, countr] + \
                                 xm[2*ii+3, countr] * (h/6) * xm[2*ii+1, countr]) * \
                                 (rhov[ii] * (vsv[ii]**2))

            # k^1 term
            snsmf[ii, countr] = snsmf[ii, countr] + \
                                (1 / (2 * (2 * np.pi * f) * U[countr])) * \
                                (xm[2*ii, countr] * (-1/2) * xm[2*ii+1, countr] + \
                                 xm[2*ii, countr] * (-1/2) * xm[2*ii+3, countr] + \
                                 xm[2*ii+1, countr] * (-1/2) * xm[2*ii, countr] + \
                                 xm[2*ii+1, countr] * (1/2) * xm[2*ii+2, countr] + \
                                 xm[2*ii+2, countr] * (1/2) * xm[2*ii+1, countr] + \
                                 xm[2*ii+2, countr] * (1/2) * xm[2*ii+3, countr] + \
                                 xm[2*ii+3, countr] * (-1/2) * xm[2*ii, countr] + \
                                 xm[2*ii+3, countr] * (1/2) * xm[2*ii+2, countr]) * \
                                 (rhov[ii] * (vsv[ii]**2))

            # k^0 term
            snsmf[ii, countr] = snsmf[ii, countr] + \
                                (vp[countr] / (2 * ((2 * np.pi * f)**2) * U[countr])) * \
                                (xm[2*ii+1, countr] * (2/h) * xm[2*ii+1, countr] + \
                                 xm[2*ii+3, countr] * (2/h) * xm[2*ii+3, countr] + \
                                 xm[2*ii+1, countr] * (-2/h) * xm[2*ii+3, countr] + \
                                 xm[2*ii+3, countr] * (-2/h) * xm[2*ii+1, countr] + \
                                 xm[2*ii, countr] * (1/h) * xm[2*ii, countr] + \
                                 xm[2*ii+2, countr] * (1/h) * xm[2*ii+2, countr] + \
                                 xm[2*ii, countr] * (-1/h) * xm[2*ii+2, countr] + \
                                 xm[2*ii+2, countr] * (-1/h) * xm[2*ii, countr]) * \
                                 (rhov[ii] * (vsv[ii]**2))

            # Thickness sensitivity
            # thickness sensitivity - omega^2 term
            snsmfh[ii, countr] = -(vp[countr] / (2 * U[countr])) * \
                                 (xm[2*ii+1, countr] * (rhov[ii]/2) * xm[2*ii+1, countr] + \
                                  xm[2*ii+3, countr] * (rhov[ii]/2) * xm[2*ii+3, countr] + \
                                  xm[2*ii+1, countr] * (rhov[ii]/2) * xm[2*ii+1, countr] + \
                                  xm[2*ii+3, countr] * (rhov[ii]/2) * xm[2*ii+3, countr]) * \
                                  h

            # k^2 term
            pmod = rhov[ii] * (vpv[ii]**2)
            smod = rhov[ii] * (vsv[ii]**2)
            snsmfh[ii, countr] = snsmfh[ii, countr] + \
                                 (1 / (2 * vp[countr] * U[countr])) * \
                                 (xm[2*ii+1, countr] * (pmod/3) * xm[2*ii+1, countr] + \
                                  xm[2*ii+3, countr] * (pmod/3) * xm[2*ii+3, countr] + \
                                  xm[2*ii+1, countr] * (pmod/6) * xm[2*ii+1, countr] + \
                                  xm[2*ii+3, countr] * (pmod/6) * xm[2*ii+3, countr] + \
                                  xm[2*ii, countr] * (smod/3) * xm[2*ii, countr] + \
                                  xm[2*ii+2, countr] * (smod/3) * xm[2*ii+2, countr] + \
                                  xm[2*ii, countr] * (smod/6) * xm[2*ii+2, countr] + \
                                  xm[2*ii+2, countr] * (smod/6) * xm[2*ii, countr]) * \
                                  h
            # k^0 term
            snsmfh[ii, countr] = snsmfh[ii, countr] + \
                                 (vp[countr] / (2 * ((2 * np.pi * f)**2) * U[countr])) * \
                                 (xm[2*ii+2, countr] * (pmod) * (-1 / (h**2)) * xm[2*ii+2, countr] + \
                                  xm[2*ii+4, countr] * (pmod) * (-1 / (h**2)) * xm[2*ii+4, countr] + \
                                  xm[2*ii+2, countr] * (-pmod) * (-1 / (h**2)) * xm[2*ii+4, countr] + \
                                  xm[2*ii+4, countr] * (-pmod) * (-1 / (h**2)) * xm[2*ii+2, countr] + \
                                  xm[2*ii+1, countr] * (smod) * (-1 / (h**2)) * xm[2*ii+1, countr] + \
                                  xm[2*ii+3, countr] * (smod) * (-1 / (h**2)) * xm[2*ii+3, countr] + \
                                  xm[2*ii+1, countr] * (-smod) * (-1 / (h**2)) * xm[2*ii+3, countr] + \
                                  xm[2*ii+3, countr] * (-smod) * (-1 / (h**2)) * xm[2*ii+1, countr]) * \
                                  h
            countr += 1
    # Special case for the bottom element
    ii = Nn-1
    h = hv[ii]
    countr = 0
    for f in fks:

        # Density sensitivity
        snsmfrho[ii, countr] = -(vp[countr] / (2 * U[countr])) * \
                               (xm[2*ii-1, countr] * (h/2) * xm[2*ii-1, countr] + \
                                xm[2*ii+0, countr] * (h/2) * xm[2*ii+0, countr]) * \
                                rhov[ii]

        # Lambda sensitivity
        # k^2 term
        snsmflam[ii, countr] = (1 / (2 * vp[countr] * U[countr])) * \
                              (xm[2*ii-1, countr] * (h/3) * xm[2*ii-1, countr]) * \
                              (rhov[ii] * ((vpv[ii]**2) - (2 * (vsv[ii]**2))))

        # k^1 term
        snsmflam[ii, countr] += (1 / (2 * (2 * np.pi * f) * U[countr])) * \
                                (xm[2*ii-1, countr] * (1/2) * xm[2*ii, countr] + \
                                 xm[2*ii, countr] * (1/2) * xm[2*ii-1, countr]) * \
                                 (rhov[ii] * ((vpv[ii]**2) - (2 * (vsv[ii]**2))))

        # k^0 term
        snsmflam[ii, countr] += (vp[countr] / (2 * ((2 * np.pi * f) ** 2) * U[countr])) * \
                                (xm[2*ii+0, countr] * (1/h) * xm[2*ii+0, countr]) * \
                                (rhov[ii] * ((vpv[ii]**2) - (2 * (vsv[ii]**2))))

        # Mu sensitivity
        # k^2 term
        snsmf[ii, countr] = (1 / (2 * vp[countr] * U[countr])) * \
                            (xm[2*ii-1, countr] * (2*h/3) * xm[2*ii-1, countr] + \
                             xm[2*ii, countr] * (h/3) * xm[2*ii, countr]) * \
                             (rhov[ii] * (vsv[ii]**2))

        # k^1 term
        snsmf[ii, countr] += (1 / (2 * (2 * np.pi * f) * U[countr])) * \
                             (xm[2*ii-1, countr] * (-1/2) * xm[2*ii, countr] + \
                              xm[2*ii, countr] * (-1/2) * xm[2*ii-1, countr]) * \
                              (rhov[ii] * (vsv[ii]**2))

        # k^0 term
        snsmf[ii, countr] += (vp[countr] / (2 * ((2 * np.pi * f) ** 2) * U[countr])) * \
                             (xm[2*ii+0, countr] * (2/h) * xm[2*ii+0, countr] + \
                              xm[2*ii-1, countr] * (1/h) * xm[2*ii-1, countr]) * \
                              (rhov[ii] * (vsv[ii]**2))

        # Thickness sensitivity
        # omega^2 term
        snsmfh[ii, countr] = -(vp[countr] / (2 * U[countr])) * \
                             (xm[2*ii-1, countr] * (rhov[ii]/2) * xm[2*ii-1, countr] + \
                              xm[2*ii+0, countr] * (rhov[ii]/2) * xm[2*ii+0, countr]) * \
                              h

        # k^2 term
        pmod = (rhov[ii] * (vpv[ii]**2))
        smod = (rhov[ii] * (vsv[ii]**2))
        snsmfh[ii, countr] += (1 / (2 * vp[countr] * U[countr])) * \
                              (xm[2*ii-1, countr] * (pmod/3) * xm[2*ii-1, countr] + \
                               xm[2*ii, countr] * (smod/3) * xm[2*ii, countr]) * \
                               h

        # k^0 term
        snsmfh[ii, countr] += (vp[countr] / (2 * ((2 * np.pi * f) ** 2) * U[countr])) * \
                              (xm[2*ii+0, countr] * (pmod) * (-1/(h**2)) * xm[2*ii+0, countr] + \
                               xm[2*ii-1, countr] * (smod) * (-1/(h**2)) * xm[2*ii-1, countr]) * \
                               h
        countr =+ 1

    # Sensitivity for Vp fixed or Vp/Vs fixed

    if pratioflag == 0:
        snsmf_vs = 2 * snsmf - 4 * np.transpose(np.transpose(snsmflam) * np.diag((vsv**2) / (vpv**2 - 2*vsv**2)))
    elif pratioflag == 1:
        snsmf_vs = 2 * (snsmf + snsmflam)
    else:
        pass
    # Make a vector of the frequencies of interest
    if vflg[0] == 0:
        vfi = [0]
    else:
        vfi = [1]

    for ii in range(1, len(vflg)):
        if vflg[ii] == 1 and vflg[ii-1] == 0:
            vfi.append(vfi[ii-1] + 2)
        elif vflg[ii] == 1 and vflg[ii-1] == 1:
            vfi.append(vfi[ii-1] + 3)
        elif vflg[ii] == 0 and vflg[ii-1] == 0:
            vfi.append(vfi[ii-1] + 1)
        else:
            vfi.append(vfi[ii-1] + 2)
    vfi=np.array(vfi)
    # Compute group kernels and change relative perturbations to absolute
    countr = 0
    Uf = np.zeros(len(fks))
    num_zeros = len(vflg) - np.count_nonzero(vflg)
    if num_zeros == 0:
        snsmf_vstot = np.zeros((Nn, len(fks))) 
        snsmf_vstotf = np.zeros((Nn, len(fks)))
        snsmf_htotf = np.zeros((Nn, len(fks)))
        snsmf_htot = np.zeros((Nn, len(fks))) 
    else :
        snsmf_vstot = np.zeros((Nn, num_zeros)) 
        snsmf_vstotf = np.zeros((Nn, num_zeros))
        snsmf_htotf = np.zeros((Nn, num_zeros))
        snsmf_htot = np.zeros((Nn, num_zeros)) 

    for f in fks_orig:
        if vflg[countr] == 0:
            Uf[countr] = vp[vfi[countr]]
            # Change phase kernel for absolute perturbation in the model and data
            # instead of relative perturbation
            snsmf_vstotf[:, countr] = np.transpose(vp[vfi[countr]] * np.transpose(snsmf_vs[:, vfi[countr]]) @ np.diag(vsv**-1))
            # Change phase kernel with respect to thickness to absolute perturbation
            snsmf_htotf[:, countr] = np.transpose(vp[vfi[countr]] * np.transpose(snsmfh[:, vfi[countr]]) @ np.diag(hv**-1))
        else:
            # The shear velocity group sensitivity kernel
            snsmf_vstot[:, countr] = snsmf_vs[:, vfi[countr]] + \
                                     ((U[vfi[countr]] / vp[vfi[countr]]) * (2 * np.pi * fks[vfi[countr]]) * \
                                     (snsmf_vs[:, vfi[countr]+1] - snsmf_vs[:, vfi[countr]-1]) / \
                                     (1 * (fks[vfi[countr]+1] - fks[vfi[countr]-1]) * 2 * np.pi))

            # For absolute perturbations
            snsmf_vstotf[:, countr] = np.transpose(vp[vfi[countr]] * np.transpose(snsmf_vstot[:, countr]) @ np.diag(vsv**-1))

            # Group velocity sensitivity kernel for changes in element thickness
            snsmf_htot[:, countr] = snsmfh[:, vfi[countr]] + \
                                    ((U[vfi[countr]] / vp[vfi[countr]]) * (2 * np.pi * fks[vfi[countr]]) * \
                                    (snsmfh[:, vfi[countr]+1] - snsmfh[:, vfi[countr]-1]) / \
                                    (1 * (fks[vfi[countr]+1] - fks[vfi[countr]-1]) * 2 * np.pi))

            # For absolute perturbations
            snsmf_htotf[:, countr] = np.transpose(U[vfi[countr]] * np.transpose(snsmf_htot[:, countr]) @ np.diag(hv**-1))
            # Decimate the group velocity for passback
            if np.isnan(U[vfi[countr]+1]) or np.isnan(U[vfi[countr]-1]):
                Uf[countr] = np.nan
            else:
                Uf[countr] = U[vfi[countr]]
        
        countr += 1



    return Uf, snsmf_vstotf, snsmf_htotf
    
def check_nans(U, U_data, fks, modn, vflg, snsmf_vstot):
    rcntr = 0
    Ur = []
    U_datar = []
    fksr = []
    fksri = []
    modnr = []
    vflgr = []
    snsmf_vstotr = []

    for ii in range(len(fks)):
        if (not np.isnan(U[ii])) and (not np.isnan(U_data[ii])):
            rcntr += 1
            Ur.append(U[ii])
            U_datar.append(U_data[ii])
            fksr.append(fks[ii])
            fksri.append(ii)
            modnr.append(modn[ii])
            vflgr.append(vflg[ii])
            snsmf_vstotr.append(snsmf_vstot[:, ii])

    Ur = np.array(Ur)
    U_datar = np.array(U_datar)
    fksr = np.array(fksr)
    fksri = np.array(fksri)
    modnr = np.array(modnr)
    vflgr = np.array(vflgr)
    snsmf_vstotr = np.column_stack(snsmf_vstotr)

    return Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr

def linvers(U_data, U, snsmf_vstot, mcmisr, dcmisr, Nn, vsv, vsg):
    # Transpose the Jacobian or kernel matrix
    G = snsmf_vstot.T
    
    # Compute the data discrepancy
    duv = (U_data - U + G @ (vsv - vsg))

    # Scale the Jacobian matrix and data
    Gs = dcmisr @ G
    duvs = dcmisr @ duv
    duvs = np.reshape(duvs, (np.size(duvs), 1))
    # Invert using lsqr
    dvs = lsqr(np.vstack((Gs, mcmisr)), np.vstack((duvs, np.zeros((Nn, 1)))), atol=1e-2, btol=1e-2)[0]#, btol=1e-2)[0]
    return dvs
    
def plot_figures(fks, U_data, U, U_guess, fksr_guess, hss, vsv_guess, vsv_update, snsmf_vstot, nupdat):
    # plot data comparisons
    U= np.trim_zeros(U, 'b')
    plt.figure()
    fsize = 16
    plt.plot(fks, U_data, "bo", linewidth=2, markersize=6)
    plt.plot(fks, U, "ko", linewidth=2, markersize=6)
    plt.plot(fksr_guess, U_guess, "ro", linewidth=2, markersize=6)
    plt.axis([1, 4, 500, 1200])
    plt.xlabel("Frequency (Hz)", fontsize=fsize, fontweight="bold")
    plt.ylabel("Velocity (m/s)", fontsize=fsize, fontweight="bold")
    plt.title("Data (blue), initial guess (red), and final update (black)")
    plt.xticks(fontsize=fsize, fontweight="bold")
    plt.yticks(fontsize=fsize, fontweight="bold")
    plt.grid(True)
    plt.tight_layout()

    # plot model comparisons
    plt.figure()
    fsize = 16
    plt.plot(vsv_guess, 0.001 * hss, "r--", linewidth=4)
    plt.plot(vsv_update[nupdat-1], 0.001 * hss, "k--", linewidth=4)
    plt.axis([700, 1500, 0, 1])
    plt.xlabel("Shear velocity (m/s)", fontsize=fsize, fontweight="bold")
    plt.ylabel("Depth (km)", fontsize=fsize, fontweight="bold")
    plt.title("Initial model (red), and final update (black)")
    plt.xticks(fontsize=fsize, fontweight="bold")
    plt.yticks(fontsize=fsize, fontweight="bold")
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.tight_layout()

####Start of driver.py ##############################################
import numpy as np
import scipy as sp
from scipy.sparse.linalg import eigs
from scipy.sparse.linalg import lsqr
from scipy.linalg import sqrtm
import matplotlib.pyplot as plt
import matplotlib
import time
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py input_value")
    else:
        target = sys.argv[1]
        # Your code to process the input_value goes here
        print(f"The input value is: {target}")
    
st = time.time()
#####################################################################
# Load input files
#####################################################################
pat='/TARGETS/Target_'+str(target)+'/'
# Load input parameters
nprm = np.loadtxt(pat+'input_params.txt')
pratioflag = nprm[0]
lsmth = nprm[1]
msigmaf = nprm[2]
nupds = int(nprm[3])
Nf = int(nprm[4])
Nn = int(nprm[5])
Nnf = int(nprm[6])
chilo = nprm[7]
chihi = nprm[8]

# Load grid
h = np.loadtxt(pat+'grid_values_solid.txt')
# Load grid in fluid
hfv = np.loadtxt(pat+'grid_values_fluid.txt')
#Load data
data = np.loadtxt(pat+'dsp_file_'+str(target)+'.dat')
# Extract the columns from the data and append to the respective lists
fks = data[:, 0]
U_data = data[:, 1]
U_data_errs = data[:, 2]
modn = data[:, 4]
vflg = data[:,3]
# Load frequencies
#fks = np.loadtxt(pat+'frequency_values.txt')
# Load velocity data
#U_data = np.loadtxt(pat+'velocity_values.txt')
# Load velocity data errors
#U_data_errs = np.loadtxt(pat+'velocity_values_errs.txt')
# Load mode numbers
#modn = np.loadtxt(pat+'mode_values.txt')
# Load velocity types
#vflg = np.loadtxt(pat+'vtype_values.txt')
# Load initial Vs model
vsv = np.loadtxt(pat+'vs_init.txt')
# Load initial Vp model
if pratioflag == 0:
    vpv = np.loadtxt(pat+'vp_init.txt')
else:
    vpvsratio = np.loadtxt(pat+'vp_init.txt')
    vpv = vpvsratio * vsv
# Load initial density model
rhov = np.loadtxt(pat+'rho_init.txt')
# Load Vp model in fluid
vpfv = np.loadtxt(pat+'vpf.txt')
# Load density model in fluid
rhofv = np.loadtxt(pat+'rhof.txt')

#####################################################################
# end input
#####################################################################

# initialize velocity update vector
vsv_update = np.zeros((nupds, Nn))

# make a vector of depths at nodes by a running sum of the grid spacings
hs = np.cumsum(h) - h[0]

# make a vector of depths at center of elements by a running sum
hss = np.cumsum(h) - h[0] + np.array(h) / 2

#####################################################################
# check whether input parameters are physically possible
#####################################################################

# density greater than zero
if np.any(rhov <= 0):
    raise ValueError('Negative density values exist in initial guess')

# shear velocity greater than zero
if np.any(vsv <= 0):
    raise ValueError('Negative shear velocity values exist in initial guess')

# Poisson's ratio between two bounds
pratio = (vpv ** 2 - 2 * (vsv ** 2)) / (2 * (vpv ** 2 - vsv ** 2))
if np.any(pratio <= -1) or np.any(pratio >= 0.5):
    raise ValueError('Impossible Poisson ratio values exist in initial guess')

# density greater than zero in fluid
if np.any(rhofv <= 0):
    raise ValueError('Negative density values exist in initial guess')
#####################################################################
# prepare for initial inversion step
#####################################################################

# Compute sensitivity kernel using initial guess
U, snsmf_vstot, snsmf_htot = raylee_sensitivity(Nn, vsv, vpv, rhov, fks, h, modn, vflg, Nnf, vpfv, rhofv, hfv, pratioflag)
# find the measurements for which both data and model are not NaN
Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr = check_nans(U, U_data, fks, modn, vflg, snsmf_vstot)
Nfr = len(fksr)
# Save the S-wave velocity guess and the resulting data
vsv_guess = vsv.copy()
U_guess = Ur.copy()
fksr_guess = fksr.copy()

# Calculate the a priori model covariance matrix and inverse square root
msigma = np.mean(U_data_errs[fksri]) * msigmaf
mcm = (msigma**2) * np.exp(-np.abs(np.tile(hs, (Nn, 1)) - np.transpose(np.tile(hs, (Nn,1)))) / lsmth)
mcmisr = sqrtm(np.linalg.inv(mcm))

# Calculate the a priori data covariance matrix and inverse square root
dcm = np.diag(U_data_errs[fksri]**2)
dcmisr = np.diag(1.0 / U_data_errs[fksri])

# Calculate the rms error of the initial guess
rmserror = np.sqrt(np.mean(((U_guess - U_datar) / 1)**2))
chisqurd = np.dot(np.dot((U_guess - U_datar), dcmisr), np.dot(dcmisr.T, (U_guess - U_datar)))
rmserror = np.atleast_1d(rmserror)
chisqurd = np.atleast_1d(chisqurd)
Nfrv = Nfr

# Check if the initial guess has chi^2 less than 1
if (chisqurd / Nfr) < chilo:
    raise ValueError('Initial model fits data to less than 1 chi-squared')
elif (chisqurd / Nfr) < chihi:
    raise ValueError('Initial model fits data within acceptable chi-squared window')
else:
    #####################################################################
    # Invert using damped least squares method of Tarantola and Valette (1982)
    #####################################################################
    # initial damped linear inversion
    
    dvs = linvers(U_datar, Ur, snsmf_vstotr, mcmisr, dcmisr, Nn, vsv, vsv_guess)
    # Add to the initial model
    vsv = dvs + vsv_guess


    if pratioflag == 1:
        vpv = vpvsratio * vsv
    else:
        pass 
    U, snsmf_vstot,_ = raylee_sensitivity(Nn, vsv, vpv, rhov, fksr, h, modnr, vflgr, Nnf, vpfv, rhofv, hfv, pratioflag)

    # Find NaNs
    U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot = check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot)
    # If the number of NaNs changed, recompute data and model covariances
    if len(fksr) != Nfr:
        Nfr = len(fksr)
        msigma = np.mean(U_data_errs[fksri]) * msigmaf
        mcm = (msigma ** 2) * np.exp(-np.abs(np.tile(hs, (Nn, 1)) - np.transpose(np.tile(hs, (Nn,1)))) / lsmth)
        mcmisr = sqrtm(np.linalg.inv(mcm))
        dcm = np.diag(U_data_errs[fksri] ** 2)
        dcmisr = np.diag(1. / U_data_errs[fksri])
    else:
        pass 
    nupdat=0
    nreds = 0
    # Compute RMS error and chi-squared
    rmserrorp = np.sqrt(np.mean(((U - U_datar) / 1) ** 2))
    chisqurdp = np.dot(np.dot((U - U_datar), dcmisr), np.dot(dcmisr.T, (U - U_datar)))
    while ((chisqurdp >= chisqurd[0] and nreds < nupds) or \
            ((chisqurdp / Nfr) < chilo and nreds < nupds)):

        nreds = nreds + 1

        # Reduce step by a factor of 2, and add it in
        dvs = dvs / 2
        vsv = vsv_update[nupdat, :] + dvs.T

        # If vpvs ratio fixed, adjust vp model
        if pratioflag == 1:
            vpv = vpvsratio * vsv
        else:
            pass
        # Call the sensitivity function to compute U
        [U, snsmf_vstot,_] = raylee_sensitivity(Nn, vsv, vpv,
                                              rhov, fksr, h, modnr, vflgr,
                                              Nnf, vpfv, rhofv, hfv, pratioflag)

        # Check for NaNs
        [U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot] = \
            check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot)

        # If the number of data points changed, adjust model and data covariances
        if len(fksr) != Nfr:
            Nfr = len(fksr)
            msigma = np.mean(U_data_errs[fksri]) * msigmaf
            mcm = (msigma ** 2) * np.exp(-np.abs(np.tile(hs, (Nn, 1)) - np.transpose(np.tile(hs, (Nn,1)))) / lsmth)
            mcmisr = sqrtm(np.linalg.inv(mcm))
            dcm = np.diag(U_data_errs[fksri] ** 2)
            dcmisr = np.diag(1. / U_data_errs[fksri])
        else:
            pass
        # The RMS of this potential update
        rmserrorp = np.sqrt(np.mean(((U - U_datar) / 1) ** 2))
        # The chi^2 of this potential update
        chisqurdp = (U - U_datar) @ dcmisr @ dcmisr.T @ (U - U_datar)
    # Shear velocity must be greater than zero
    if np.sum(vsv <= 0) > 0:
        raise ValueError('Negative shear velocity values encountered in inversion')
    else:
        pass

    # Poisson's ratio between two bounds
    pratio = (vpv**2 - 2 * vsv**2) / (2 * (vpv**2 - vsv**2))
    if np.sum(pratio <= -1) > 0 or np.sum(pratio >= 0.5) > 0:
        raise ValueError('Impossible Poisson ratio values encountered in inversion')
    else:
        pass
    nupdat=0
    vsv_update[nupdat, :] = vsv
    # The rms of this update
    rmserror = np.append(rmserror,rmserrorp)
    # The chi^2 of this update
    chisqurd = np.append(chisqurd,chisqurdp)
# Full modeling
[U, snsmf_vstot,_] = raylee_sensitivity(Nn, vsv, vpv, rhov, fks, h, modn, vflg, Nnf, vpfv, rhofv, hfv, pratioflag)

# Check for NaNs
[Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr] = check_nans(U, U_data, fks, modn, vflg, snsmf_vstot)
# If number of NaNs changed, recompute data and model covariances
if len(fksr) != Nfr:
    Nfr = len(fksr)
    msigma = np.mean(U_data_errs[fksri]) * msigmaf
    mcm = (msigma**2) * np.exp(-np.abs(np.tile(hs, (Nn, 1)) - np.transpose(np.tile(hs, (Nn,1)))) / lsmth)
    mcmisr = sqrtm(np.linalg.inv(mcm))
    dcm = np.diag(U_data_errs[fksri]**2)
    dcmisr = np.diag(1.0 / U_data_errs[fksri])
else:
    pass
# Compute RMS and chi-squared
rmserrorp = np.sqrt(np.mean(((Ur - U_datar) / 1) ** 2))
chisqurdp = np.dot(np.dot((Ur - U_datar), dcmisr), np.dot(dcmisr.T, (Ur - U_datar)))
# The rms of this update
rmserror = np.append(rmserror,rmserrorp)
# The chi^2 of this update
chisqurd = np.append(chisqurd,chisqurdp)
Nfrv =  np.append(Nfrv, Nfr)

# Now an iterative loop, updating the initial guess
# While the stopping criterion and the maximum
# allowed number of iterations have not been met, continue updating
while chisqurdp / Nfr > chihi and nupdat < nupds:
    # Perform updates here
    

    # Invert again as in Tarantola and Valette (1982)
    # Linear inverse
    dvs = linvers(U_datar, Ur, snsmf_vstotr, mcmisr, dcmisr, Nn, vsv, vsv_guess)

    # Add to the initial model
    vsv = dvs + vsv_guess
    # If fixed vpvs ratio, adjust vp model
    if pratioflag == 1:
        vpv = vpvsratio * vsv
    # Call the sensitivity function to model
    U, snsmf_vstot,_ = raylee_sensitivity(Nn, vsv, vpv, rhov, fksr, h, modnr, vflgr, Nnf, vpfv, rhofv, hfv, pratioflag)

    # Check for NaNs
    U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot = check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot)
    # If the number of data changed, recompute data and model covariances
    if len(fksr) != Nfr:
        Nfr = len(fksr)
        msigma = np.mean(U_data_errs[fksri]) * msigmaf
        mcm = (msigma ** 2) * np.exp(-np.abs(np.tile(hs, (Nn, 1)) - np.transpose(np.tile(hs, (Nn,1)))) / lsmth)
        mcmisr = sqrtm(np.linalg.inv(mcm))
        dcm = np.diag(U_data_errs[fksri] ** 2)
        dcmisr = np.diag(1. / U_data_errs[fksri])
    else:
        pass
    # Compute RMS and chi of this potential update
    rmserrorp = np.sqrt(np.mean(((U - U_datar) / 1) ** 2))
    chisqurdp = np.dot(np.dot((U - U_datar), dcmisr), np.dot(dcmisr.T, (U - U_datar)))

    # A reduced line search if chi^2 of update is not lower
    nreds = 0
    # The gradient - difference between the current update and previous
    dvs = (vsv - vsv_update[nupdat, :]).T
    while ((chisqurdp >= 1.01 * chisqurd[nupdat] and nreds < nupds) or \
            ((chisqurdp / Nfr) < chilo and nreds < nupds)):

        nreds = nreds + 1

        # Reduce step by a factor of 2, and add it in
        dvs = dvs / 2
        vsv = vsv_update[nupdat, :] + dvs.T

        # If vpvs ratio fixed, adjust vp model
        if pratioflag == 1:
            vpv = vpvsratio * vsv
        else:
            pass
        # Call the sensitivity function to compute U
        [U, snsmf_vstot,_] = raylee_sensitivity(Nn, vsv, vpv,
                                              rhov, fksr, h, modnr, vflgr,
                                              Nnf, vpfv, rhofv, hfv, pratioflag)

        # Check for NaNs
        [U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot] = \
            check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot)

        # If the number of data points changed, adjust model and data covariances
        if len(fksr) != Nfr:
            Nfr = len(fksr)
            msigma = np.mean(U_data_errs[fksri]) * msigmaf
            mcm = (msigma ** 2) * np.exp(-np.abs(np.tile(hs, (Nn, 1)) - np.transpose(np.tile(hs, (Nn,1)))) / lsmth)
            mcmisr = np.linalg.sqrtm(np.linalg.inv(mcm))
            dcm = np.diag(U_data_errs(fksri) ** 2)
            dcmisr = np.diag(1. / U_data_errs[fksri])
        else:
            pass
        # The RMS of this potential update
        rmserrorp = np.sqrt(np.mean(((U - U_datar) / 1) ** 2))
        # The chi^2 of this potential update
        chisqurdp = (U - U_datar) @ dcmisr @ dcmisr.T @ (U - U_datar)
    
    # Shear velocity must be greater than zero
    if np.sum(vsv <= 0) > 0:
        raise ValueError('Negative shear velocity values encountered in inversion')
    else:
        pass

    # Poisson's ratio between two bounds
    pratio = (vpv**2 - 2 * vsv**2) / (2 * (vpv**2 - vsv**2))
    if np.sum(pratio <= -1) > 0 or np.sum(pratio >= 0.5) > 0:
        raise ValueError('Impossible Poisson ratio values encountered in inversion')
    else:
        pass
    # The next updated model, print number of update to screen
    vsv_update[nupdat, :] = vsv
    # The rms of this update
    rmserror = np.append(rmserror,rmserrorp)
    # The chi^2 of this update
    chisqurd = np.append(chisqurd,chisqurdp)
    
    # Now full modeling
    U, snsmf_vstot,_ = raylee_sensitivity(Nn, vsv, vpv, rhov, fks, h, modn, vflg, Nnf, vpfv, rhofv, hfv, pratioflag)
    # Check for NaNs
    Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr = check_nans(U, U_data, fks, modn, vflg, snsmf_vstot)

    # If the number of data points changed, recompute data and model covariances
    if len(fksr) != Nfr:
        Nfr = len(fksr)
        msigma = np.mean(U_data_errs[fksri]) * msigmaf
        mcm = (msigma**2) * np.exp(-np.abs(np.tile(hs, (Nn, 1)) - np.transpose(np.tile(hs, (Nn,1)))) / lsmth)
        mcmisr = sqrtm(np.linalg.inv(mcm))
        dcm = np.diag(U_data_errs[fksri]**2)
        dcmisr = np.diag(1. / U_data_errs[fksri])
    else:
        pass
    # Compute RMS and chi^2
    rmserrorp = np.sqrt(np.mean(((Ur - U_datar) / 1) ** 2))
    chisqurdp = np.dot(np.dot(Ur - U_datar, dcmisr), dcmisr.T).dot(Ur - U_datar)

    # The rms of this update
    rmserror = np.append(rmserror,rmserrorp)
    # The chi^2 of this update
    chisqurd = np.append(chisqurd,chisqurdp)
    Nfrv =  np.append(Nfrv, Nfr)
    
    # Increment the update counter
    nupdat += 1
    
#plot_figures(fks, U_data, U, U_guess, fksr_guess, hss, vsv_guess, vsv_update, snsmf_vstot, nupdat)
combined_vs = np.column_stack((hss / 1000, vsv_update[nupdat-1]))
np.savetxt(pat+"final_vs.txt", combined_vs, fmt='%.8f')
U= np.trim_zeros(U, 'b')
fw_md = np.column_stack((fks,U))
np.savetxt(pat+"fw_group.txt", fw_md, fmt='%.8f')

print(f"{Nfr} of {Nf - np.isnan(U_data).sum()} measurements used")
et = time.time()
# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
if chisqurd[nupdat + 1] / Nfr > chihi:
    "WARNING: Inversion did not converge to stopping criterion and underfitted data. Increase number of updates."
else:
    pass

if chisqurd[nupdat + 1] / Nfr < chilo:
    "WARNING: Inversion did not converge to stopping criterion and overfitted data. Increase number of reduction steps."
else:
    pass
