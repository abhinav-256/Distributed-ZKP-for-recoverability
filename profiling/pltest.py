from phe import paillier
public_key, private_key = paillier.generate_paillier_keypair()
for i in range(50):
    encrypted_number=public_key.encrypt(3)
    decrypted_number=private_key.decrypt(encrypted_number)

