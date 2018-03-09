︠4fad2ef6-cc8a-48c4-a942-b61553007cfb︠
""" Sagemath worksheet for Advanced Crypto Ex2
Questions 1b, and 2
Guy Coop - 1447634 - gtc434 """

import hashlib


# provided values
p = 19557005198709002903944860293031210393112763219086296641657077101075696217230957575932867920574816565856404011416908876733063974554042055901449481509845537165833627083609136426386625960849573203936690977282752217830467255076903813033768202255463477317279507253420233052258146393806024179130344501759742626253484354878829682461413002143954786668386100202257237258007726886368228024987052808444032800050116408995869561195208576647050961019393809940446490709061546503439488751331252440431940049716566765873615456614911458658341732120931319694200689362617296630728389003180670326151653001322354906363854478585334821235883

q = Integer((p-1)/2)

g1 = 4;
g2 = 3
x = ZZ.random_element(q);
g = g1;
pk = power_mod(g,x,p);

# mathematical methods
def div_mod(numerator, denominator, mod):
    return (numerator * pow(denominator, -1, mod)) % mod

# elgamal methods
def elgamal_enc(p, q, g, pk, r, m):
    R = pow(g, r, p)
    S = (pow(pk, r, p) * pow(g, m, p)) % p
    return (R, S)

# other methods
def get_hashed(values):
    hash_string = ""
    for value in values:
        hash_string += str(value)
    hash_value = hashlib.sha256(hash_string).hexdigest()
    return int(hash_value, 16)

# question 1b implement DisjProof = (DisjProve, DisjVerify)
def DisjProve(p, q, g, pk, R, S, r, m):
    if m == 0:
        u_0 = ZZ.random_element(q); U_0 = pow(g, u_0, p); V_0 = pow(pk, u_0, p) # real proof
        ch_1 = ZZ.random_element(q); rep_1 = ZZ.random_element(q) # fake proof
        U_1 = div_mod(pow(g, rep_1, p), pow(R, ch_1, p), p) # fake proof
        V_1 = div_mod(pow(pk, rep_1, p), pow(S * pow(g, -1, p), ch_1, p), p) # fake proof
        ch = get_hashed([R, S, U_0, V_0, U_1, V_1]) # real proof
        ch_0 = ch - ch_1; rep_0 = (u_0 + (ch_0 * r)) # disjunction
    else:
        ch_0 = ZZ.random_element(q); rep_0 = ZZ.random_element(q) # fake proof
        U_0 = div_mod(pow(g, rep_0, p), pow(R, ch_0, p), p) # fake proof
        V_0 = div_mod(pow(pk, rep_0, p), pow(S, ch_0, p), p) # fake proof
        u_1 = ZZ.random_element(q); U_1 = pow(g, u_1, p); V_1 = pow(pk, u_1, p) # real proof
        ch = get_hashed([R, S, U_0, V_0, U_1, V_1]) # real proof
        ch_1 = ch - ch_0; rep_1 = (u_1 + (ch_1 * r))  # disjunction
    return [ch_0, ch_1, rep_0, rep_1]


def DisjVerify(p, q, g, pk, R, S, ch_0, ch_1, rep_0, rep_1):
    hashed_value = get_hashed([
            R, S,
            div_mod(pow(g, rep_0, p), pow(R, ch_0, p), p),
            div_mod(pow(pk, rep_0, p), pow(S, ch_0, p), p),
            div_mod(pow(g, rep_1, p), pow(R, ch_1, p), p),
            div_mod(pow(pk, rep_1, p), pow(S * pow(g, -1, p), ch_1, p), p)
        ])
    return ch_0 + ch_1 == hashed_value

# test DisjProve and DisjVerify
m = 0
r = ZZ.random_element(q)
[R, S] = elgamal_enc(p, q, g, pk, r, m)
[ch_0, ch_1, rep_0, rep_1] = DisjProve(p, q, g, pk, R, S, r, m)
print( "m = 0 gives:{}".format(DisjVerify(p, q, g, pk, R, S, ch_0, ch_1, rep_0, rep_1)) )

m = 1
r = ZZ.random_element(q)
[R, S] = elgamal_enc(p, q, g, pk, r, m)
[ch_0, ch_1, rep_0, rep_1] = DisjProve(p, q, g, pk, R, S, r, m)
print( "m = 1 gives:{}".format(DisjVerify(p, q, g, pk, R, S, ch_0, ch_1, rep_0, rep_1)) )

m = 2
r = ZZ.random_element(q)
[R, S] = elgamal_enc(p, q, g, pk, r, m)
[ch_0, ch_1, rep_0, rep_1] = DisjProve(p, q, g, pk, R, S, r, m)
print( "m = 2 gives:{}".format(DisjVerify(p, q, g, pk, R, S, ch_0, ch_1, rep_0, rep_1)) )



# provided methods
def eqdlprove(p,q,g1,g2,y1,y2,x):
  r = ZZ.random_element(q);
  com1 = power_mod(g1,r,p);
  com2 = power_mod(g2,r,p);
  hash_value = hashlib.sha256(str(y1) + str(y2) + str(com1) + str(com2)).hexdigest();
  ch = Integer(int(hash_value,16),q);
  s = Integer(Mod(r + x * ch,q));
  return [ch,s]

# implement question 2
def setup(g, p, q):
    sk = ZZ.random_element(q)
    pk = pow(g, sk, p)
    return pk, sk


def vote(p, q, pk, v):
    r = ZZ.random_element(q)
    R, S = elgamal_enc(p, q, g, pk, r, v)
    [ch_0, ch_1, rep_0, rep_1] = DisjProve(p, q, g, pk, R, S, r, m)
    return [R, S], [ch_0, ch_1, rep_0, rep_1]


def verify(p, q, g, b):
    pk = b[0]
    [R, S] = b[1]
    [ch_0, ch_1, rep_0, rep_1] = b[2]
    return DisjVerify(p, q, g, pk, R, S, ch_0, ch_1, rep_0, rep_1)


def tally(p, q, g, BB, sk)
    # in the description, i cant understand where pk camr from
    pk = pow(g, sk, p)
    R_sum = 1
    S_sum = 1
    for ballot in BB:
        if not verify(p, q, g, ballot):
            return False, None #return result_valid = false, and no count
        [R, S] = ballot[1]
        R_sum = (R_sum * R) % p
        S_sum = (S_sum * S) % p
    result = S_sum * pow(pow(R_sum, sk, p), -1, p)
    proof = eqdlprove(p, q, g, pk, pow(result, -1, p), sk)
    return result, proof




︡fe43ed33-d2c6-4649-b33d-9df3ed36b7f4︡{"stdout":"' Sagemath worksheet for Advanced Crypto Ex2\\nQuestions 1b, and 2\\nGuy Coop - 1447634 - gtc434 '\n"}︡{"stdout":"m = 0 gives:True\n"}︡{"stdout":"m = 1 gives:True\n"}︡{"stdout":"m = 2 gives:False\n"}︡{"stdout":"False\n"}︡{"done":true}︡
︠8c898c9d-3f8c-436d-94dd-ebc3eb408236︠









