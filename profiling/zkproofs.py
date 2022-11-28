from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,G2,GT,pair
import random

group = PairingGroup('SS1024')
f,g,h = [group.random(G1) for i in range(3)]

#### Utils ####

def sz(base64_input):
    b_padded = base64_input.split(str.encode(":"))[1]
    pad_size = b_padded.count(str.encode("="))
    b_len_without_pad = len(b_padded)-4
    byte_len = (b_len_without_pad *3)/4 +(3-pad_size)-1
    bit_len = byte_len * 8
    return bit_len

#### Commitments ####

def commit(m, _r):
	return (g**m) * (h**_r)

def open(c, m, _r):
	return c == (g**m) * (h**_r)

def test_commitment():
	m = 3
	_r = group.random(ZR)
	c = commit(m, _r)
	cserial = group.serialize(c)
	cprime = group.deserialize(cserial)
	return open(cprime, m, _r)


#### BLS signatures ####

def blskeygen():
	_sk = group.random(ZR)
	pk = g ** _sk
	return _sk, pk

def blssign(m, _sk):
	mhash = group.hash(m, type=G1)
	sigma = mhash ** _sk
	return sigma

def blsverify(sigma, m, pk):
	mhash = group.hash(m, type=G1)
	return pair(mhash, pk) == pair(sigma, g)

def test_bls():
	_sk, pk = blskeygen()
	m = b'Hello'
	sigma = blssign(m, _sk)
	return blsverify(sigma, m, pk)


#### BB signatures ####

def bbkeygen():
	_sk = group.random(ZR)
	pk = g ** _sk
	return _sk, pk

def bbsign(m, _sk):
	sigma = g ** (1 / (m + _sk))
	return sigma

def bbverify(sigma, m, pk):
	return pair(sigma, pk * (g ** m)) == pair(g, g)

def test_bb():
	_sk, pk = bbkeygen()
	m = 5
	sigma = bbsign(m, _sk)
	return bbverify(sigma, m, pk)


#### BBS+ signatures ####

def bbspluskeygen():
	_sk = group.random(ZR)
	pk = f ** _sk
	return _sk, pk

def bbsplussign(m, _sk):
	c, r = group.random(ZR, 2)
	A = (f * (g ** m) * (h ** r))** (1 / (c + _sk))
	return (A, c, r)

def bbsplusverify(sigma, m, pk):
	A, c, r = sigma
	return pair(A, pk * (f ** c)) == pair(f * (g ** m) * (h ** r), f)

def bbsplussign_commitment(C, _sk):
	c, rdash = group.random(ZR, 2)
	A = (f * (h ** rdash) * C)** (1 / (c + _sk))
	return (A, c, rdash)

def bbsplussign_randomise(sigma, r):
	A, c, rdash = sigma
	rdoubledash = rdash + r
	return (A, c, rdoubledash)

def test_bbsplus():
	_sk, pk = bbspluskeygen()
	m = 5
	sigma = bbsplussign(m, _sk)
	return bbsplusverify(sigma, m, pk)

def test_bbsplus_commitment():
	_sk, pk = bbspluskeygen()
	m = 5
	r = group.random(ZR)
	C = commit(m, r)
	sigma = bbsplussign_commitment(C, _sk)
	sigmasdoubledash = bbsplussign_randomise(sigma, r)
	return bbsplusverify(sigmasdoubledash, m, pk)


#### PoK1: PoK of a signature on a committed message ####
#### PK{(rho, r, v): C = g^rho h^r and V = g^{v/(x+rho)}}

def pok1commit(V):
	_s, _t, _m = group.random(ZR, 3)
	a = (pair(V, g) ** (-_s)) * (pair(g, g) ** _t)
	D = (g ** _s) * (h ** _m)
	return _s, _t, _m, a, D

def pok1chal():
	return group.random(ZR)

def pok1resp(_s, _t, _m, _rho, _r, _v, c):
	zrho = _s - (_rho * c)
	zv = _t - (_v * c)
	zr = _m - (_r * c)
	return zrho, zv, zr

def pok1verify(C, V, a, D, zrho, zv, zr, y, c):
	cond1 = (D == (C ** c) * (h ** zr) * (g ** zrho))
	cond2 = (a == (pair(V, y) ** c) * (pair(V, g) ** (-zrho)) * (pair(g, g) ** (zv)))
	return cond1 and cond2 

def pok1(C, V, _rho, _r, _v, y):
	_s, _t, _m, a, D = pok1commit(V)
	c = pok1chal()
	zrho, zv, zr = pok1resp(_s, _t, _m, _rho, _r, _v, c)
	return pok1verify(C, V, a, D, zrho, zv, zr, y, c)

def pok1nizkproof(C, V, _rho, _r, _v, y):
	stmt = (g, h, C, V)
	_s, _t, _m, a, D = pok1commit(V)
	c = group.hash((stmt, (a, D)), type=ZR)
	zrho, zv, zr = pok1resp(_s, _t, _m, _rho, _r, _v, c)
	return c, (zrho, zv, zr)

def pok1nizkverify(C, V, y, pf):
	stmt = (g, h, C, V)
	c, (zrho, zv, zr) = pf
	verif = ((C ** c) * (h ** zr) * (g ** zrho), (pair(V, y) ** c) * (pair(V, g) ** (-zrho)) * (pair(g, g) ** (zv)))
	return c == group.hash((stmt, verif), type=ZR)


#### Set membership proof (C commits some rho in set phi) ####

def genphi(n):
	return [group.random(ZR) for i in range(n)]

def gencomms(phi):
	n = len(phi)
	_shift = random.randint(0,n-1)
	_rands = []
	comms = []
	for i in range(n):
		j = (i + _shift) % n
		_rand = group.random(ZR)
		comm = commit(phi[j], _rand)
		_rands.append(_rand)
		comms.append(comm)
	return _shift, _rands, comms

def verfsigs(phi):
	_sk, pk = bbkeygen()
	sigs = []
	for m in phi:
		sig = bbsign(m, _sk)
		sigs.append(sig)
	return pk, sigs

def test_set_membership():
	n = 5
	phi = genphi(n)
	_shift, _rands, comms = gencomms(phi)

	# Verifier decides to check i^th commitment in comms
	i = random.randint(0,n-1)
	C = comms[i]

	# Verifier sends fresh public key and signatures on all elements of phi
	group.InitBenchmark()
	group.StartBenchmark(["RealTime"])
	y, sigs = verfsigs(phi)
	group.EndBenchmark()
	print(" - Generating %s verifier signatures (sec):" % n, group.GetBenchmark("RealTime"))

	# Prover identifies the randomness and the signature corresponding to the i^th commitment
	_r = _rands[i]
	_j = (i + _shift) % n
	_sig = sigs[_j]
	_rho = phi[_j]
	
	# Prover sends a blinded signature on the shifted message
	_v = group.random(ZR)
	V = _sig ** _v

	# Proof of knowledge that a signature on a committed message is known
	group.StartBenchmark(["RealTime"])
	status1 = pok1(C, V, _rho, _r, _v, y)
	group.EndBenchmark()
	print(" - Proving for a single element (sec):", group.GetBenchmark("RealTime"))

	# NIZK version of the above proof
	group.StartBenchmark(["RealTime"])
	pf = pok1nizkproof(C, V, _rho, _r, _v, y)
	group.EndBenchmark()
	print(" - Generating NIZK proof for a single element (sec):", group.GetBenchmark("RealTime"))
	c, (zrho, zv, zr) = pf
	print(" - Size of NIZK proof for a single element (bits):", sz(group.serialize(c)) + sz(group.serialize(zrho)) + sz(group.serialize(zv)) + sz(group.serialize(zr)))
	group.StartBenchmark(["RealTime"])
	status2 = pok1nizkverify(C, V, y, pf)
	group.EndBenchmark()
	print(" - Verifying NIZK proof for a single element (sec):", group.GetBenchmark("RealTime"))

	return status1 and status2

#### POK2: PoK of a representation of a group element ####
#### PK{(rho,r): C = g^rho h^r}

def pok2commit():
	_r1, _r2 = group.random(ZR, 2)
	a = (g ** _r1) * (h ** _r2)
	return _r1, _r2, a

def pok2resp(_r1, _r2, _rho, _r, c):
	z1 = _r1 - (_rho * c)
	z2 = _r2 - (_r * c)
	return z1, z2

def pok2nizkproof(C, _rho, _r):
	stmt = (g, h, C)
	_r1, _r2, a = pok2commit()
	c = group.hash((stmt, a), type=ZR)
	z1, z2 = pok2resp(_r1, _r2, _rho, _r, c)
	return c, (z1, z2)

def pok2nizkverify(C, pf):
	stmt = (g, h, C)
	c, (z1, z2) = pf
	verif = (C ** c) * (g ** z1) * (h ** z2)
	return c == group.hash((stmt, verif), type=ZR)

def pok2nizkproofs(comms, phi, _shift, _rands):
	n = len(comms)
	pfs = []
	for i,comm in enumerate(comms):
		_j = (i + _shift) % n
		_rho = phi[_j]
		_r = _rands[i]
		pf = pok2nizkproof(comm, _rho, _r)
		pfs.append(pf)
	return pfs

def pok2nizkverifies(comms, pfs):
	status = True
	for i,comm in enumerate(comms):
		pf = pfs[i]
		status = status and pok2nizkverify(comm, pf)
	return status

def test_pok2():
	_rho, _r = group.random(ZR, 2)
	C = commit(_rho, _r)
	pf = pok2nizkproof(C, _rho, _r)
	return pok2nizkverify(C, pf)

#### PoK3: PoK of a BBS+ signature on a given message ####
#### PK{(s1, s2, c, r, d1, d2): B1 = g^s1 h^s2 and B1^c = g^d1 h^d2 and
####       e(B2, y)e(f, f)^(-1) = e(B2, f)^{-c}e(g, y)^s2 e(g, f)^d2 e(g, f)^rho e(h, f)^r}

def pok3nizkproof(B1, B2, y, rho, _s1, _s2, _c, _r, _d1, _d2):
	iden = g ** 0
	invB1 = B1 ** (-1)
	eg1 = pair(B2, y)
	eg2 = pair(f, f)
	eg3 = pair(B2, f)
	inveg3 = eg3 ** (-1)
	eg4 = pair(g, y)
	eg5 = pair(g, f)
	eg6 = pair(h, f)
	egstar = eg1 * (eg2 ** (-1)) * (eg5 ** (-rho))

	# Now proving PK{(s1, s2, c, r, d1, d2): 
	#                  B1     = g^s1 h^s2             AND 
	#                  iden   = invB1^c g^d1 h^d2     AND 
	#                  egstar = inveg3^c eg4^s2 eg5^d2 eg6^r }

	stmt = (g, h, invB1, inveg3, eg4, eg5, eg6, B1, iden, egstar)
	
	# Commit (each ri corresponds to the respective known discrete log)
	_r1, _r2, _r3, _r4, _r5, _r6 = group.random(ZR, 6)
	com1 = (g ** _r1) * (h ** _r2)
	com2 = (invB1 ** _r3) * (g ** _r5) * (h ** _r6)
	com3 = (inveg3 ** _r3) * (eg4 ** _r2) * (eg5 ** _r6) * (eg6 ** _r4)

	# Challenge
	c = group.hash((stmt, (com1, com2, com3)), type=ZR)

	# Response
	z1 = _r1 - c * _s1
	z2 = _r2 - c * _s2
	z3 = _r3 - c * _c 
	z4 = _r4 - c * _r
	z5 = _r5 - c * _d1
	z6 = _r6 - c * _d2

	return c, (z1, z2, z3, z4, z5, z6)

def pok3nizkverify(B1, B2, y, rho, pf):
	iden = g ** 0
	invB1 = B1 ** (-1)
	eg1 = pair(B2, y)
	eg2 = pair(f, f)
	eg3 = pair(B2, f)
	inveg3 = eg3 ** (-1)
	eg4 = pair(g, y)
	eg5 = pair(g, f)
	eg6 = pair(h, f)
	egstar = eg1 * (eg2 ** (-1)) * (eg5 ** (-rho))

	# Now verifying PK{(s1, s2, c, r, d1, d2): 
	#                  B1     = g^s1 h^s2             AND 
	#                  iden   = invB1^c g^d1 h^d2     AND 
	#                  egstar = inveg3^c eg4^s2 eg5^d2 eg6^r }

	stmt = (g, h, invB1, inveg3, eg4, eg5, eg6, B1, iden, egstar)

	c, (z1, z2, z3, z4, z5, z6) = pf
	
	v1 = (B1 ** c) * (g ** z1) * (h ** z2)
	v2 = (iden ** c) * (invB1 ** z3) * (g ** z5) * (h ** z6)
	v3 = (egstar ** c) * (inveg3 ** z3) * (eg4 ** z2) * (eg5 ** z6) * (eg6 ** z4)

	return c == group.hash((stmt, (v1, v2, v3)), type=ZR)

def test_pok3():
	rho = group.random(ZR)
	_x, y = bbspluskeygen()
	_A, _c, _r = bbsplussign(rho, _x)
	_s1, _s2 = group.random(ZR, 2)
	_d1 = _c * _s1
	_d2 = _c * _s2
	B1 = commit(_s1, _s2)
	B2 = _A * (g ** _s2)
	pf = pok3nizkproof(B1, B2, y, rho, _s1, _s2, _c, _r, _d1, _d2)
	return pok3nizkverify(B1, B2, y, rho, pf)

#### Reverse set membership proof (rho is committed by some C in set phi') ####
#### PK{(r): g^rho h^r \in Phi'}

def verfsigs_commitment(comms):
	_sk, pk = bbspluskeygen()
	sigs = []
	for comm in comms:
		sig = bbsplussign_commitment(comm, _sk)
		sigs.append(sig)
	return pk, sigs

def test_reverse_set_membership():
	n = 5
	phi = genphi(n)
	_shift, _rands, comms = gencomms(phi)

	# Prover proves knowledge of the committed values
	group.InitBenchmark()
	group.StartBenchmark(["RealTime"])
	pfs = pok2nizkproofs(comms, phi, _shift, _rands)
	status = pok2nizkverifies(comms, pfs)
	group.EndBenchmark()
	print(" - Proving knowledge of %s committed values (sec):" % n, group.GetBenchmark("RealTime"))

	# Verifier sends fresh public key and signatures on all elements of comms
	group.InitBenchmark()
	group.StartBenchmark(["RealTime"])
	y, sigs = verfsigs_commitment(comms)
	group.EndBenchmark()
	print(" - Generating %s verifier signatures (sec):" % n, group.GetBenchmark("RealTime"))

	# Verifier decides to check j^th element in phi
	j = random.randint(0,n-1)
	rho = phi[j]

	# Prover identifies the index for the commitment committing rho
	_i = (j - _shift) % n
	_sig = sigs[_i]
	_rand = _rands[_i]
	
	# NIZK proof generation
	group.StartBenchmark(["RealTime"])
	_sigdoubledash = bbsplussign_randomise(_sig, _rand)
	_A, _c, _r = _sigdoubledash
	_s1, _s2 = group.random(ZR, 2)
	_d1 = _c * _s1
	_d2 = _c * _s2
	B1 = commit(_s1, _s2)
	B2 = _A * (g ** _s2)
	pf = pok3nizkproof(B1, B2, y, rho, _s1, _s2, _c, _r, _d1, _d2)
	group.EndBenchmark()
	print(" - Generating NIZK proof for a single element (sec):", group.GetBenchmark("RealTime"))
	c, (z1, z2, z3, z4, z5, z6) = pf
	print(" - Size of NIZK proof for a single element (bits):", sz(group.serialize(c)) + sz(group.serialize(z1)) + sz(group.serialize(z2)) + sz(group.serialize(z3)) + sz(group.serialize(z4)) + sz(group.serialize(z5)) + sz(group.serialize(z6)))
	
	# NIZK proof verification
	group.StartBenchmark(["RealTime"])
	status = pok3nizkverify(B1, B2, y, rho, pf)
	group.EndBenchmark()
	print(" - Verifying NIZK proof for a single element (sec):", group.GetBenchmark("RealTime"))

	return status


if __name__ == "__main__":
	print("Commitments ok?:", test_commitment())
	print("BLS signatures ok?:", test_bls())
	print("BB signatures ok?:", test_bb())
	print("BBS+ signatures ok?:", test_bbsplus())
	print("BBS+ signatures on committed messages ok?:", test_bbsplus_commitment())
	print("Set membership ok?:", test_set_membership())
	print("Knowledge of representation of group element ok?:", test_pok2())
	print("Knowledge of BBS+ signature ok?:", test_pok3())
	print("Reverse set membership ok?:", test_reverse_set_membership())
