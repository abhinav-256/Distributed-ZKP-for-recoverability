from elgamal import elgamal_keygen_par, elgamal_encrypt, \
    elgamal_encrypt_with_randomness, elgamal_decrypt, elgamal_reencrypt,elgamal_reencrypt_with_randomness, \
    elgamal_mult,elgamal_div, elgamal_exp
from perm import permute_par, gen_rand_perm_par
from charm.toolbox.pairinggroup import ZR
from globals import group, g1
import random

def check(_authpk,randomness,cp1,cp2):
    status_permute = ( elgamal_div(cp1,cp2) == elgamal_encrypt_with_randomness(_authpk, 1, randomness))
    return status_permute 

def enc_sigs_rev_permute(sigs,myn,_authpk,_authsk,alpha,pi,re_pi):    
    # Encrypt sigs under elgamal encryption
    enc_sigs = [0]*myn
    elgamal_randomness = [0]*myn
    for i in range(len(sigs)):
        elgamal_randomness[i] = group.random(ZR)
        enc_sigs[i] = elgamal_encrypt_with_randomness(_authpk, sigs[i], elgamal_randomness[i])

    # Re-encrypt and reverse-permute signatures
    store_enc_sigs=[]
    store_elgamal_randomness=[]
    store_enc_sigs.append(enc_sigs)

    for a in reversed(range(alpha)):
        reenc_sigs =[0]*len(enc_sigs)
        reenc_randomness =[0]*len(enc_sigs)
        for i in range(len(enc_sigs)):
            reenc_randomness[i] = group.random(ZR)
            reenc_sigs[i] = elgamal_reencrypt_with_randomness(_authpk, enc_sigs[i], reenc_randomness[i] )
        enc_sigs   = permute_par(reenc_sigs, re_pi[a])
        enc_randomness = permute_par(reenc_randomness, re_pi[a])
        store_enc_sigs.append(enc_sigs)
        store_elgamal_randomness.append(enc_randomness)

    store_enc_sigs = list(reversed(store_enc_sigs))
    store_elgamal_randomness = list(reversed(store_elgamal_randomness))
    status_permute =True 
    for a in range(1, alpha-1, 2):
        audit_a = [random.randint(0,1) for _ in range(myn)]
        for i in range(myn):
            if audit_a[i] == 0: # RIGHT
                ii = pi[a][i]
                status_permute = status_permute and check(_authpk,store_elgamal_randomness[a][i],store_enc_sigs[a][i], store_enc_sigs[a+1][ii])
            else: # LEFT
                ii = re_pi[a-1][i]
                status_permute = status_permute and check(_authpk,store_elgamal_randomness[a-1][ii],store_enc_sigs[a-1][ii], store_enc_sigs[a][i])

    print("Status Elgamal Permute:",(status_permute))
    return enc_sigs