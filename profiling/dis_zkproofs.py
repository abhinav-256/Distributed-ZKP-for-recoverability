import os
import random
import sys

from charm.toolbox.pairinggroup import G1, ZR

from globals import g1, group, pai_group, beta, q, n, alpha
from parutils import nprocs, init_strios, print_strios, close_strios, pprint, fullname
from pedersen import commit
from bbsig import bbkeygen_par, bbsign, bbverify, bbbatchverify
from bbsplussig import bbspluskeygen_par, bbsplussign, bbsplusverify, bbsplussign_commitment, bbsplussign_obtain, bbsplusquasibatchverify
from paillier import pai_keygen_par, pai_encrypt, pai_decrypt, pai_reencrypt, pai_add
from elgamal import elgamal_keygen_par, elgamal_encrypt, elgamal_decrypt, elgamal_reencrypt, elgamal_mult, elgamal_exp 
from secretsharing import gen_beaver_triples
from perm import permute_par, gen_rand_perm_par
from pok import dpk0nizkproofs, dpk0nizkverifs, dpk1nizkproofs, dpk1nizkverifs, dpk2nizkproofs, dpk2nizkverifs
from misc import statusstr, fmt, timed, timedperEA, timer, timerperEA, retval, sz

random.seed()

try:
    os.environ['verbosity'] = sys.argv[3]
except Exception:
    os.environ['verbosity'] = '0'

def ballotgen_par(n, alpha):
    """ Generate commitments, vote shares and randomness shares 

    Note: For benchmarking purposes, only generating one (commitment, vote 
    share, randomness share) tuple per ballot, representing the voter's 
    choice, as opposed to m such tuples. """

    myn = n//nprocs

    vote_shares = [None]*alpha
    randomness_shares = [None]*alpha
    comm_shares = [None]*alpha
    for a in range(alpha):
        vote_shares[a] = [group.random(ZR) for _ in range(myn)]
        randomness_shares[a] = [group.random(ZR) for _ in range(myn)]
        comm_shares[a] = \
            [commit(vote_shares[a][i], randomness_shares[a][i]) for i in range(myn)]

    comms = [g1 ** 0] * myn
    for i in range(myn):
        for a in range(alpha):
            comms[i] = comms[i] * comm_shares[a][i]
    
    return comms, vote_shares, randomness_shares

def genperms_par(n, alpha):
    """ Generate secret permutations (and reverse permutations) for each authority. """

    pi = []
    re_pi =[]
    for a in range(alpha):
        pi_a, re_pi_a = gen_rand_perm_par(n)
        pi.append(pi_a)
        re_pi.append(re_pi_a)
    return pi, re_pi

@timedperEA
def process_votes(vote_shares, pi, alpha):
    """ Process votes using the stored vote shares and permutations at each authority. """

    myn = len(vote_shares[0])
    _auth_paisk, auth_paipk = pai_keygen_par()

    # Generate encrypted vote shares
    enc_vote_shares = []
    for a in range(alpha):
        vote_shares_a = vote_shares[a]
        enc_vote_shares_a = []
        for i in range(myn):
            enc_vote_share_a = pai_encrypt(auth_paipk, vote_shares_a[i]) 
            enc_vote_shares_a.append(enc_vote_share_a)
        enc_vote_shares.append(enc_vote_shares_a)

    # Combine encrypted vote shares
    enc_votes = [pai_encrypt(auth_paipk, 0)] * myn
    for a in range(alpha):
        enc_vote_shares_a = enc_vote_shares[a]
        for i in range(myn):
            enc_votes[i] = pai_add(auth_paipk, enc_votes[i], enc_vote_shares_a[i])

    # Re-encrypt and permute 
    for a in range(alpha):
        reenc_votes = [pai_reencrypt(auth_paipk, enc_vote) for enc_vote in enc_votes]
        enc_votes = permute_par(reenc_votes, pi[a])

    # Add integer blinding
    for a in range(alpha):
        for i in range(myn):
            b_a = q*random.randint(0, beta)
            enc_b_a  = pai_encrypt(auth_paipk, b_a)
            enc_votes[i] = pai_add(auth_paipk, enc_votes[i], enc_b_a)

    # Decryption
    votes_out = []
    for i in range(myn):
        pai_vote = pai_decrypt(auth_paipk, _auth_paisk, enc_votes[i])
        vote = group.init(ZR, int(pai_vote))
        votes_out.append(vote)

    return votes_out

@timed
def get_verfsigs(votes_out):
    """ Get verifier signatures on published votes. """

    _verfsk, verfpk = bbkeygen_par()

    return verfpk, [bbsign(vote_out, _verfsk) for vote_out in votes_out]

@timedperEA
@retval
def check_verfsigs(sigs, votes_out, verfpk):
    """ Check whether verifier signatures on published votes are correct. """

    return bbbatchverify(sigs, votes_out, verfpk)

@timedperEA
def get_blsigs(sigs, re_pi, alpha):
    """ Obtain blinded signatures on the published vote outputs, but ordered by the same
    order as the published commitments. """

    myn = len(sigs)

    _authsk, authpk = elgamal_keygen_par()

    # Encrypt sigs under elgamal encryption
    with timerperEA("encrypt sigs", verbosity=1):
        enc_sigs = [0]*myn
        for i in range(len(sigs)):
            enc_sigs[i] = elgamal_encrypt(authpk, sigs[i])

    # Re-encrypt and reverse-permute signatures
    with timerperEA("re-encrypt and reverse-permute signatures", verbosity=1):
        for a in reversed(range(alpha)):
            reenc_sigs = [elgamal_reencrypt(authpk, enc_sig) for enc_sig in enc_sigs]
            enc_sigs   = permute_par(reenc_sigs, re_pi[a])

    # Generate blinding shares
    with timerperEA("getting blinded shares", verbosity=1):
        blshares = []
        for a in range(alpha):
            blshares_a = [group.random(ZR) for _ in range(n)]
            blshares.append(blshares_a)
    
    # Generate encrypted blinded signatures
    with timerperEA("generate encrypted blinded signatures", verbosity=1): 
        enc_blsigs = [elgamal_encrypt(authpk, 1)] * myn
        for a in range(alpha):
            enc_blsigs_a = [0]*myn
            for i in range(myn):
                enc_blsigs_a[i] = elgamal_exp(enc_sigs[i], blshares[a][i])
                enc_blsigs[i] = elgamal_mult(enc_blsigs[i], enc_blsigs_a[i])

    # Decryption of blinded signatures
    with timerperEA("decrypt blinded signatures", verbosity=1): 
        blsigs = [elgamal_decrypt(_authsk, enc_blsig) for enc_blsig in enc_blsigs]

    return blsigs, blshares

@timed
def get_verfsigs_rev(comms):
    _verfsk, verfpk = bbspluskeygen_par()
    return verfpk, [bbsplussign_commitment(comm, _verfsk) for comm in comms]

@timedperEA
@retval
def check_verfsigs_rev(sigs, comms, verfpk):
    """ Check whether verifier quasi-signatures on published commitments are correct. """

    return bbsplusquasibatchverify(sigs, comms, verfpk)

@timedperEA
def get_blsigs_rev(sigs, randomness_shares, pi, alpha):
    """ Obtain blinded signatures on the published vote outputs for the 
    reverse set membership proof. """

    myn = len(sigs)

    _auth_elgsk, auth_elgpk = elgamal_keygen_par()
    _auth_paisk, auth_paipk = pai_keygen_par()

    # Encrypt components of sigs
    with timerperEA("encrypt sigs", verbosity=1):
        enc_sigs_S    = [g1 ** 0] * myn
        enc_sigs_c    = [0] * myn
        enc_sigs_rbar = [0] * myn
        for i in range(len(sigs)):
            (S, c, rbar) = sigs[i]
            enc_sigs_S[i] = elgamal_encrypt(auth_elgpk, S)
            enc_sigs_c[i] = pai_encrypt(auth_paipk, c)
            enc_sigs_rbar[i] = pai_encrypt(auth_paipk, rbar)   
        
    # Encrypt all randomness
    with timerperEA("encrypt randomness (paillier)", verbosity=1):
        enc_randomness =[]
        for a in range(alpha):
                randomness_shares_a = randomness_shares[a]
                enc_randomness_a = [pai_encrypt(auth_paipk, randomness_share_a) for randomness_share_a in randomness_shares_a]
                enc_randomness.append(enc_randomness_a)

    # Homomorphically obtain the r component of the BBS+ signature
    with timerperEA("homomorphically obtain r component of BBS+ sig", verbosity=1):
        enc_sigs_r = enc_sigs_rbar
        for a in range(alpha):
            enc_randomness_a = enc_randomness[a]
            for i in range(myn):
                enc_sigs_r[i] = pai_add(auth_paipk, enc_sigs_r[i], enc_randomness_a[i])

    # Re-encrypt and permute signatures
    with timerperEA("re-encrypt and permute", verbosity=1):
        for a in range(alpha):
            reenc_sigs_S = [elgamal_reencrypt(auth_elgpk, enc_sig_S) for enc_sig_S in enc_sigs_S]
            reenc_sigs_c = [pai_reencrypt(auth_paipk, enc_sig_c) for enc_sig_c in enc_sigs_c]
            reenc_sigs_r = [pai_reencrypt(auth_paipk, enc_sig_r) for enc_sig_r in enc_sigs_r]
            enc_sigs_S = permute_par(reenc_sigs_S, pi[a])
            enc_sigs_c = permute_par(reenc_sigs_c, pi[a])
            enc_sigs_r = permute_par(reenc_sigs_r, pi[a])

    # Generate blinding shares
    with timerperEA("generate blinding shares", verbosity=1):
        blshares_S = []
        blshares_c = []
        blshares_r = []
        for a in range(alpha):
            blshares_S_a = [group.random(ZR) for _ in range(myn)]
            blshares_c_a = [group.init(ZR, random.randint(0, beta*q)) for _ in range(myn)]
            blshares_r_a = [group.init(ZR, random.randint(0, beta*q)) for _ in range(myn)]
            blshares_S.append(blshares_S_a)
            blshares_c.append(blshares_c_a)
            blshares_r.append(blshares_r_a)

    # Generate encrypted blinded signatures
    with timerperEA("generate encrypted blinded signatures", verbosity=1):
        enc_blsigs_S = enc_sigs_S
        enc_blsigs_c = enc_sigs_c
        enc_blsigs_r = enc_sigs_r
        for a in range(alpha):
            enc_blshares_S_a = [0] * myn
            enc_blshares_c_a = [0] * myn
            enc_blshares_r_a = [0] * myn
            for i in range(myn):
                enc_blshares_S_a[i] = elgamal_encrypt(auth_elgpk, g1**blshares_S[a][i])
                enc_blshares_c_a[i] = pai_encrypt(auth_paipk, blshares_c[a][i])
                enc_blshares_r_a[i] = pai_encrypt(auth_paipk, blshares_r[a][i])
                
                enc_blsigs_S[i] = elgamal_mult(enc_blsigs_S[i], enc_blshares_S_a[i])
                enc_blsigs_c[i] = pai_add(auth_paipk, enc_blsigs_c[i], enc_blshares_c_a[i])
                enc_blsigs_r[i] = pai_add(auth_paipk, enc_blsigs_r[i], enc_blshares_r_a[i])

    # Decryption of blinded signatures
    with timerperEA("decryption of blinded signatures", verbosity=1):
        blsigs_S = [elgamal_decrypt(_auth_elgsk, enc_blsig_S) for enc_blsig_S in enc_blsigs_S]
        blsigs_c = [group.init(ZR, int(pai_decrypt(auth_paipk, _auth_paisk, enc_blsig_c))) for enc_blsig_c in enc_blsigs_c]
        blsigs_r = [group.init(ZR, int(pai_decrypt(auth_paipk, _auth_paisk, enc_blsig_r))) for enc_blsig_r in enc_blsigs_r]

    return (blsigs_S, blsigs_c, blsigs_r), (blshares_S, blshares_c, blshares_r) 

@timed
def main(n, alpha):
    """ Main benchmarking code for processing and proving forward and reverse set membership 
    for n votes with alpha authorities. """

    pai_sk, pai_pk = pai_keygen_par()
    elg_sk, elg_pk = elgamal_keygen_par()

    comms, vote_shares, randomness_shares = ballotgen_par(n, alpha)
    pi, re_pi = genperms_par(n, alpha)
    votes_out = process_votes(vote_shares, pi, alpha)
    

    ####### Forward set membership proof ########

    # Get blinded signatures and blinded shares
    verfpk, sigs = get_verfsigs(votes_out)
    status_verfsigs_fwd = check_verfsigs(sigs, votes_out, verfpk)
    blsigs, blshares = get_blsigs(sigs, re_pi, alpha)

    # Proofs
    dpk0_pfs = dpk0nizkproofs(comms, blsigs, vote_shares, randomness_shares, blshares, verfpk, alpha)

    # Verification
    status_fwd = dpk0nizkverifs(comms, blsigs, dpk0_pfs, verfpk)

    ####### Reverse set membership proof ########

    # Generate beaver triples (offline preprocessing step for multiplicative 
    # secret sharing used in DPK2)
    gen_beaver_triples(2*len(votes_out), alpha)

    # Prover proves knowledge of the committed values
    dpk1_pfs = dpk1nizkproofs(comms, vote_shares, randomness_shares, alpha)
    status_comm = dpk1nizkverifs(comms, dpk1_pfs)

    # Verifier sends fresh public key and BBS+ signatures
    verfpk, sigs_rev = get_verfsigs_rev(comms)
    status_verfsigs_rev = check_verfsigs_rev(sigs_rev, comms, verfpk)
    blsigs_rev, blshares_rev = get_blsigs_rev(sigs_rev, randomness_shares, pi, alpha)

    # Proofs
    blsigs_S, blsigs_c, blsigs_r = blsigs_rev
    blshares_S, blshares_c, blshares_r = blshares_rev
    dpk2_pfs = dpk2nizkproofs(votes_out, blsigs_S, blsigs_c, blsigs_r, blshares_S, blshares_c, blshares_r, verfpk, alpha)

    # Correctness of proof
    status_rev = dpk2nizkverifs(votes_out, blsigs_S, blsigs_c, blsigs_r, dpk2_pfs, verfpk)

    # Sizes
    pprint("size of one group element in G1:", sz(group.serialize(group.random(G1))), " bytes")
    pprint("size of one group element in Zq:", sz(group.serialize(group.random(ZR))), " bytes")
    pprint("size of one BB signature:", sz(group.serialize(sigs[0])), "bytes")
    pprint("size of one BBS+ signature:", sz(group.serialize(sigs_rev[0][0])) + sz(group.serialize(sigs_rev[0][1])) + sz(group.serialize(sigs_rev[0][2])), "bytes")
    pprint("size of one DPK0 proof:", sum([sz(group.serialize(pf_item)) for pf_item in dpk0_pfs[0]]), "bytes")
    pprint("size of one DPK1 proof:", sum([sz(group.serialize(pf_item)) for pf_item in dpk1_pfs[0]]), "bytes")
    pprint("size of one DPK2 proof:", sum([sz(group.serialize(pf_item)) for pf_item in dpk2_pfs[0]]), "bytes")


if __name__ == "__main__":
    try:
        init_strios()
        main(n, alpha)
    finally:
        print_strios()
        close_strios()
