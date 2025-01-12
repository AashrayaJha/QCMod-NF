R = Integers(27)
A = matrix([[R(0), 26], [4, 6]])
B = matrix([[R(10), 1], [25, 26]])
H = [matrix.identity(2, R)]
for M in H:
     if not A*M in H:
           H.append(A*M)
     if not B*M in H:
         H.append(B*M)
assert(len(H) == 1944)
P = matrix([[R(0), 1], [-1, 0]])
P*A.transpose()*P^(-1) in H
#returns True
P*B.transpose()*P^(-1) in H
#returns True