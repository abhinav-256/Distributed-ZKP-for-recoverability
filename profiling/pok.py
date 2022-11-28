from charm.toolbox.pairinggroup import ZR, pair

from globals import group, g1, h1, f1, f2, fT, eg1f2, inveg1f2, inveh1f2, idenT
from secretsharing import sharerand, sharemult
from misc import timed, timedperEA, retval


#### DPK 0: PK{(v,r,bl): C = g1^v h1^r AND e(blsig^(1/bl), yf2^v) = e(g1, f2)} ########

def dpk0nizkproof(C, blsig, _v, _r, _bl, verfpk, alpha):
    """ PoK of a signature on a committed message, i.e.:
     PK{(v, r, bl): C = g^v h^r and blsig^{1/bl} = g^{1/(x+v)}} """

    # Commit
    R = eg1f2 ** 0
    D = g1 ** 0
    _s, _t, _m = [0]*alpha, [0]*alpha, [0]*alpha
    for a in range(alpha):
        _s[a], _t[a], _m[a] = group.random(ZR, 3)
        R_a = (pair(blsig, f2 ** (-_s[a]))) * (eg1f2 ** _t[a])
        D_a = (g1 ** _s[a]) * (h1 ** _m[a])        
        R = R * R_a
        D = D * D_a

    # Challenge

    c = group.hash((g1, h1, f2, C, blsig) + (R, D), type=ZR)

    # Response
    zv = group.init(ZR, 0)
    zbl = group.init(ZR, 0)
    zr = group.init(ZR, 0)
    for a in range(alpha):
        zv_a = _s[a] - (_v[a] * c)
        zbl_a = _t[a] - (_bl[a] * c)
        zr_a = _m[a] - (_r[a] * c)
        zv = zv + zv_a
        zbl = zbl + zbl_a
        zr = zr + zr_a

    return c, zv, zbl, zr

def dpk0nizkverify(C, blsig, verfpk, pf):
    stmt = (g1, h1, f2, C, blsig)

    c, zv, zbl, zr = pf
    verif = ((pair(blsig, (verfpk**c) * (f2**(-zv))) * (eg1f2 ** (zbl))),
             (C ** c) * (h1 ** zr) * (g1 ** zv))
    return c == group.hash((g1, h1, f2, C, blsig) + verif, type=ZR)

@timedperEA
def dpk0nizkproofs(comms, blsigs, vote_shares, randomness_shares, blind_shares, verfpk, alpha):
    pfs = []
    for i in range(len(comms)):
        _v, _r, _bl = [],[],[]
        for a in range(alpha):
            _v.append(vote_shares[a][i])
            _r.append(randomness_shares[a][i])
            _bl.append(blind_shares[a][i])

        # NIZK proof
        pf = dpk0nizkproof(comms[i], blsigs[i], _v, _r, _bl, verfpk, alpha)
        pfs.append(pf)
    return pfs

@timed
@retval
def dpk0nizkverifs(comms, blsigs, pfs, verfpk):
    status = True
    for i in range(len(pfs)):
        status = status and dpk0nizkverify(comms[i], blsigs[i], verfpk, pfs[i])
    return status

#### DPK 1: PK{(v,r): C = g1^v h1^r} ############################################

def dpk1nizkproof(C, _v, _r, alpha):
    """ PoK of the opening of the commitment, i.e.:
     PK{(rho, r): C = g^rho h^r} """
    stmt = (g1, h1, C)

    # Commit
    s = g1 ** 0
    _r1, _r2 = [0]*alpha,[0]*alpha
    for a in range(alpha):
        _r1[a], _r2[a] = group.random(ZR, 2)
        s_a = (g1 ** _r1[a]) * (h1 ** _r2[a])
        s = s * s_a

    # Challenge
    c = group.hash(stmt+ (s,), type=ZR)
   
    # Response
    z1 = group.init(ZR, 0)
    z2 = group.init(ZR, 0)
    for a in range(alpha):
        z1_a = _r1[a] - (_v[a] * c)
        z2_a = _r2[a] - (_r[a] * c)
        z1 =z1+z1_a
        z2 =z2+z2_a

    return c, z1, z2

def dpk1nizkverify(C, pf):
    stmt = (g1, h1, C)
    c, z1, z2 = pf
    verif = (C ** c) * (g1 ** z1) * (h1 ** z2)
    
    return c == group.hash(stmt +(verif,), type=ZR)

@timedperEA
def dpk1nizkproofs(comms, vote_shares, randomness_shares, alpha):
    pfs = []
    for i in range(len(comms)):
        _v, _r = [],[]
        for a in range(alpha):
            _v.append(vote_shares[a][i])
            _r.append(randomness_shares[a][i])

        pf = dpk1nizkproof(comms[i], _v, _r, alpha)
        pfs.append(pf)
    return pfs

@timed
@retval
def dpk1nizkverifs(comms, pfs):
    status = True
    for i in range(len(pfs)):
        status = status and dpk1nizkverify(comms[i], pfs[i])
    return status

#### DPK 2: PK{(bS,bC,br,delta0,delta1,delta2):               #######################################
####                  eta   = g1^bc g2^bS g3^br g4^delta1 AND #######################################
####                  z     = g4^bS g5^delta0             AND #######################################
####                  idenT = z^-bc g4^delta1 g5^delta2 }     #######################################

def dpk2nizkproof(v, Stilde, ctilde, rtilde, _bS, _bc, _br, verfpk, alpha):
    eta  = pair(Stilde, verfpk * (f2 ** ctilde)) / pair(f1 * (g1 ** v) * (h1 ** rtilde), f2)
    gen1 = pair(Stilde, f2)
    gen2 = pair(g1, verfpk * (f2 ** ctilde))
    gen3 = inveh1f2
    gen4 = inveg1f2
    gen5 = fT

    # Generate shares of delta0, delta1, delta2. 
    _delta0 = sharerand(alpha)
    _delta1 = sharemult(_bS, _bc, alpha)
    _delta2 = sharemult(_delta0, _bc, alpha)

    z = eg1f2 ** 0
    for a in range(alpha):
        z_a = (gen4 ** _bS[a]) * (gen5 ** _delta0[a])
        z = z * z_a 

    stmt = (f1, g1, h1, f2, Stilde, ctilde, rtilde, 
            idenT, eta, gen1, gen2, gen3, gen4, gen5, z) 

    # Commit
    com1 = eg1f2 ** 0
    com2 = eg1f2 ** 0
    com3 = eg1f2 ** 0
    _rbS, _rbc, _rbr, _rdelta0, _rdelta1, _rdelta2 = [0]*alpha, [0]*alpha, [0]*alpha, [0]*alpha, [0]*alpha, [0]*alpha
    for a in range(alpha):
        _rbS[a], _rbc[a], _rbr[a], _rdelta0[a], _rdelta1[a], _rdelta2[a] = group.random(ZR, 6)
        com1_a = (gen1 ** _rbc[a]) * (gen2 ** _rbS[a]) * (gen3 ** _rbr[a]) * (gen4 ** _rdelta1[a])
        com2_a = (gen4 ** _rbS[a]) * (gen5 ** _rdelta0[a])
        com3_a = (z ** (-_rbc[a])) * (gen4 ** _rdelta1[a]) * (gen5 * _rdelta2[a])
        com1 = com1 * com1_a
        com2 = com2 * com2_a
        com3 = com3 * com3_a

    # Challenge
    c = group.hash((stmt + (com1, com2, com3)), type=ZR)

    # Response
    zbS     = group.init(ZR, 0)
    zbc     = group.init(ZR, 0)
    zbr     = group.init(ZR, 0)
    zdelta0 = group.init(ZR, 0)
    zdelta1 = group.init(ZR, 0)
    zdelta2 = group.init(ZR, 0)
    for a in range(alpha):
        zbS_a     = _rbS[a]     - c * _bS[a]
        zbc_a     = _rbc[a]     - c * _bc[a]
        zbr_a     = _rbr[a]     - c * _br[a]
        zdelta0_a = _rdelta0[a] - c * _delta0[a]
        zdelta1_a = _rdelta1[a] - c * _delta1[a]
        zdelta2_a = _rdelta2[a] - c * _delta2[a]
        
        zbS       = zbS + zbS_a
        zbc       = zbc + zbc_a
        zbr       = zbr + zbr_a
        zdelta0   = zdelta0 + zdelta0_a
        zdelta1   = zdelta1 + zdelta1_a
        zdelta2   = zdelta2 + zdelta2_a

    return z, c, zbS, zbc, zbr, zdelta0, zdelta1, zdelta2
    

def dpk2nizkverify(v, Stilde, ctilde, rtilde, pf, verfpk):
    eta = pair(Stilde, verfpk * (f2 ** ctilde)) / pair(f1 * (g1 ** v) * (h1 ** rtilde), f2)
    gen1 = pair(Stilde, f2)
    gen2 = pair(g1, verfpk * (f2 ** ctilde))
    gen3 = inveh1f2
    gen4 = inveg1f2
    gen5 = fT

    z, c, zbS, zbc, zbr, zdelta0, zdelta1, zdelta2 = pf
    stmt = (f1, g1, h1, f2, Stilde, ctilde, rtilde, 
            idenT, eta, gen1, gen2, gen3, gen4, gen5, z)


    v1 = (eta ** c) * (gen1 ** zbc) * (gen2 ** zbS) * (gen3 ** zbr) * (gen4 ** zdelta1)
    v2 = (z ** c) * (gen4 ** zbS) * (gen5 ** zdelta0)
    v3 = (idenT ** c) * (z ** (-zbc)) * (gen4 ** zdelta1) * (gen5 * zdelta2)

    return c == group.hash((stmt + (v1, v2, v3)), type=ZR)

@timedperEA
def dpk2nizkproofs(votes, blsigs_S, blsigs_c, blsigs_r, blind_shares_S, blind_shares_c, blind_shares_r, verfpk, alpha):
    pfs = []
    for j in range(len(votes)):
        vote = votes[j]

        _bS, _bc, _br = [],[],[]
        for a in range(alpha):
            _bS.append(blind_shares_S[a][j])
            _bc.append(blind_shares_c[a][j])
            _br.append(blind_shares_r[a][j])

        pf = dpk2nizkproof(votes[j], blsigs_S[j], blsigs_c[j], blsigs_r[j],  _bS, _bc, _br, verfpk, alpha)
        pfs.append(pf)
    return pfs

@timed
@retval
def dpk2nizkverifs(votes, blsigs_S, blsigs_c, blsigs_r, pfs, verfpk):
    status = True
    for j in range(len(votes)):
        status = status and dpk2nizkverify(votes[j], blsigs_S[j], blsigs_c[j], blsigs_r[j], pfs[j], verfpk)
    return status
