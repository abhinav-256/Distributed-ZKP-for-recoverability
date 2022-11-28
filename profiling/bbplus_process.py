from elgamal import elgamal_keygen_par, elgamal_encrypt, \
    elgamal_encrypt_with_randomness, elgamal_decrypt, elgamal_reencrypt,elgamal_reencrypt_with_randomness, \
    elgamal_mult,elgamal_div, elgamal_exp
from paillier import pai_keygen_par, pai_encrypt, pai_decrypt, pai_reencrypt, pai_add, pai_encrypt_with_randomness, pai_reencrypt_with_randomness, pai_div
from perm import permute_par, gen_rand_perm_par
from charm.toolbox.pairinggroup import ZR
from globals import g1, group, pai_group, beta, q
import random

def check(auth_elgpk,auth_paipk,elg_randomness,elg_cp1,elg_cp2,pai1_randomness,pai1_cp1,pai1_cp2,pai2_randomness,pai2_cp1,pai2_cp2):
    #check for corretness of elgamal, pailler 1 and paillier 2
    bbplus_status_permute = ( elgamal_div(elg_cp1,elg_cp2) == elgamal_encrypt_with_randomness(auth_elgpk, 1, elg_randomness)) and (pai_div(auth_paipk,pai1_cp1,pai1_cp2) == pai_encrypt_with_randomness(auth_paipk, 0,pai1_randomness)) and (pai_div(auth_paipk,pai2_cp1,pai2_cp2) == pai_encrypt_with_randomness(auth_paipk, 0,pai2_randomness))
    return bbplus_status_permute 

# Re-encrypt and reverse-permute signatures
def reenc_repermute_bbplus(myn, alpha, auth_paipk,auth_elgpk,enc_sigs_S,enc_sigs_c,enc_sigs_r,pi,re_pi): 
    store_enc_sigs_S=[]
    store_enc_sigs_c=[]
    store_enc_sigs_r=[]
    store_enc_sigs_S.append(enc_sigs_S)
    store_enc_sigs_c.append(enc_sigs_c)
    store_enc_sigs_r.append(enc_sigs_r)
    store_elgamal_randomness_1=[] 
    store_paillier_randomness_1=[]
    store_paillier_randomness_2=[]
    for a in range(alpha):
        reenc_sigs_S = [0]*len(enc_sigs_S)
        reenc_sigs_c = [0]*len(enc_sigs_S)
        reenc_sigs_r = [0]*len(enc_sigs_S)
        reenc_elgamal_randomness1=[0]*len(enc_sigs_S)
        reenc_paillier_randomness1=[0]*len(enc_sigs_S)
        reenc_paillier_randomness2=[0]*len(enc_sigs_S)
        for i in range(len(enc_sigs_S)):
            reenc_elgamal_randomness1[i] = group.random(ZR)
            reenc_paillier_randomness1[i] = pai_group.random(auth_paipk[0])
            reenc_paillier_randomness2[i] = pai_group.random(auth_paipk[0])
            reenc_sigs_S[i] = elgamal_reencrypt_with_randomness(auth_elgpk, enc_sigs_S[i],reenc_elgamal_randomness1[i])
            reenc_sigs_c[i] = pai_reencrypt_with_randomness(auth_paipk, enc_sigs_c[i],reenc_paillier_randomness1[i])
            reenc_sigs_r[i] = pai_reencrypt_with_randomness(auth_paipk, enc_sigs_r[i],reenc_paillier_randomness2[i])
        enc_sigs_S = permute_par(reenc_sigs_S, pi[a])
        enc_sigs_c = permute_par(reenc_sigs_c, pi[a])
        enc_sigs_r = permute_par(reenc_sigs_r, pi[a])
        enc_elgamal_randomess1 = permute_par(reenc_elgamal_randomness1, pi[a])
        enc_paillier_randomess1 = permute_par(reenc_paillier_randomness1, pi[a])
        enc_paillier_randomess2 = permute_par(reenc_paillier_randomness2, pi[a])
        store_enc_sigs_S.append(enc_sigs_S)
        store_enc_sigs_c.append(enc_sigs_c)
        store_enc_sigs_r.append(enc_sigs_r)
        store_elgamal_randomness_1.append(enc_elgamal_randomess1)
        store_paillier_randomness_1.append(enc_paillier_randomess1)
        store_paillier_randomness_2.append(enc_paillier_randomess2)

    bbplus_status_permute =True 
    for a in range(1, alpha-1, 2):
        audit_a = [random.randint(0,1) for _ in range(myn)]
        for i in range(myn):
            if audit_a[i] == 0: # RIGHT
                ii = pi[a][i]
                bbplus_status_permute = bbplus_status_permute and check(auth_elgpk,auth_paipk,store_elgamal_randomness_1[a][ii],store_enc_sigs_S[a+1][ii], store_enc_sigs_S[a][i], store_paillier_randomness_1[a][ii],store_enc_sigs_c[a+1][ii], store_enc_sigs_c[a][i], store_paillier_randomness_2[a][ii],store_enc_sigs_r[a+1][ii], store_enc_sigs_r[a][i])
                                        
            else: # LEFT
                ii = re_pi[a-1][i]
                bbplus_status_permute = bbplus_status_permute and check(auth_elgpk,auth_paipk,store_elgamal_randomness_1[a-1][i],store_enc_sigs_S[a][i], store_enc_sigs_S[a-1][ii],store_paillier_randomness_1[a-1][i],store_enc_sigs_c[a][i], store_enc_sigs_c[a-1][ii],store_paillier_randomness_2[a-1][i],store_enc_sigs_r[a][i], store_enc_sigs_r[a-1][ii])
    print("reverse", bbplus_status_permute)
    return enc_sigs_S, enc_sigs_c,enc_sigs_r