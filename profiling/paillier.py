# Reference --https://jhuisi.github.io/charm/charm/schemes/pkenc/pkenc_paillier99.html


from charm.toolbox.integergroup import lcm, integer, toInt

from globals import pai_group, pai_p, pai_q, pai_n 
from parutils import rank,serialize_bcast


def L(pai_u, pai_n):
    """ computes L(u) => ((u - 1) / n) """
    pai_U = integer(int(pai_u) - 1)
    if int(pai_U) == 0:
        return integer(0, pai_n)
    return pai_U / pai_n

def pai_keygen():
    lam = lcm(pai_p - 1, pai_q - 1)
    pai_n2 = pai_n ** 2
    pai_g = pai_group.random(pai_n2)
    pai_u = (L(((pai_g % pai_n2) ** lam), pai_n) % pai_n) ** -1
    pai_pk, pai_sk = [[pai_n, pai_g, pai_n2], [lam, pai_u]]
    return (pai_sk, pai_pk)

def pai_add(pai_pk, cipher1, cipher2):
    return (cipher1*cipher2) % pai_pk[2]

def pai_div(pai_pk,cipher1, cipher2):
    return (cipher1/cipher2) % pai_pk[2]

def pai_encrypt(pai_pk, pai_m):
    pai_g, pai_n, pai_n2 = pai_pk[1], pai_pk[0], pai_pk[2]
    pai_r = pai_group.random(pai_pk[0])
    pai_c = ((pai_g % pai_n2) ** int(pai_m)) * ((pai_r % pai_n2) ** pai_n)
    return pai_c

def pai_encrypt_with_randomness(pai_pk, pai_m,pai_r):
    pai_g, pai_n, pai_n2 = pai_pk[1], pai_pk[0], pai_pk[2]
    pai_c = ((pai_g % pai_n2) ** int(pai_m)) * ((pai_r % pai_n2) ** pai_n)
    return pai_c

def pai_decrypt(pai_pk, pai_sk, pai_c):
    pai_n, pai_n2 = pai_pk[0], pai_pk[2]
    pai_m = ((L(pai_c ** pai_sk[0], pai_n) %pai_n) * pai_sk[1]) % pai_n
    return toInt(pai_m)

def pai_reencrypt_with_randomness(pai_pk, cipher,r):
    return pai_add(pai_pk, cipher, pai_encrypt_with_randomness(pai_pk, 0,r))

def pai_reencrypt(pai_pk, cipher):
    return pai_add(pai_pk, cipher, pai_encrypt(pai_pk, 0))

def pai_keygen_par():
    if rank == 0:
        pai_sk, pai_pk = pai_keygen()
    else:
        pai_sk, pai_pk = None, None
    mypai_sk, mypai_pk = serialize_bcast(pai_sk, pai_pk)
    return mypai_sk, mypai_pk

import sys

if __name__ == "__main__":
    sk, pk = pai_keygen()
    for i in range(int(sys.argv[1])):
        c = pai_encrypt(pk, 3)
        m = pai_decrypt(pk, sk, c)
    print("m=", m)
