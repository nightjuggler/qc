#!/usr/bin/python
#
# qc.py - Quantum Computing Library
#
#         Various functions for mathematically simulating the quantum circuit model
#         of quantum computation including examples of superdense coding and quantum
#         teleportation.
#
# Version 1.0 - February 13-16, 2016 - Pius Fischer
#
import math
import random

def validState(V, minSize=2):
	vlen = len(V)
	assert vlen & (vlen - 1) == 0
	assert vlen >= minSize
	return vlen

def validMatrix(U, size=None):
	numRows = len(U)

	if size is None:
		assert numRows & (numRows - 1) == 0
	else:
		assert numRows == size

	for row in U:
		assert len(row) == numRows

	return numRows

def multiplyMatrixByScalar(scalar, U):
	validMatrix(U)

	return [[scalar * element for element in row] for row in U]

def multiplyMatrixByMatrix(U1, U2):
	validMatrix(U2, validMatrix(U1))

	U2 = zip(*U2)
	return [[sum([e1 * e2 for e1, e2 in zip(row1, row2)]) for row2 in U2] for row1 in U1]

def changeState(U, V):
	# Input:
	#   U is a unitary matrix (n x n)
	#   V is a state vector (n x 1)
	# Output:
	#   V which is overwritten by the product of U and V

	validMatrix(U, validState(V))

	newState = [sum([ue * ve for ue, ve in zip(row, V)]) for row in U]

	for i, element in enumerate(newState):
		V[i] = element

	return V

def combineStates(V1, V2):
	# Input:
	#   V1 is a state vector (m x 1)
	#   V2 is a state vector (n x 1)
	# Output:
	#   A state vector (mn x 1) that is the Kronecker product of V1 and V2

	validState(V1)
	validState(V2)

	return [e1 * e2 for e1 in V1 for e2 in V2]

def combineTransforms(U1, U2):
	# Input:
	#   U1 is a unitary matrix (m x m)
	#   U2 is a unitary matrix (n x n)
	# Output:
	#   A unitary matrix (mn x mn) that is the Kronecker product of U1 and U2

	validMatrix(U1)
	validMatrix(U2)

	return [[e1 * e2 for e1 in row1 for e2 in row2] for row1 in U1 for row2 in U2]

IdentityMatrix = ((1, 0), (0, 1))

def PhaseShiftGate(phi):
	return ((1, 0), (0, math.cos(phi) + math.sin(phi) * 1J))

def ControlledGate(U):
	validMatrix(U, 2)

	return ((1, 0, 0, 0),
		(0, 1, 0, 0),
		(0, 0, U[0][0], U[0][1]),
		(0, 0, U[1][0], U[1][1]))

HadamardGate = multiplyMatrixByScalar(1/math.sqrt(2), ((1, 1), (1, -1)))

XGate = ((0, 1), (1, 0))
YGate = ((0, -1J), (1J, 0))
ZGate = ((1, 0), (0, -1))

ControlledNotGate = ControlledGate(XGate)

def qubitLength(V):
	vlen = validState(V)
	n = 0
	while vlen > 1:
		n += 1
		vlen >>= 1
	return n

def printQubit(V):
	n = qubitLength(V)

	print '['
	for i, e in enumerate(V):
		print ('  {:0' + str(n) + 'b} -> {: }  p={}').format(i, e, abs(e) * abs(e))
	print ']'

def reorderState(V, i):
	# Input:
	#   V is a state vector for n qubits (n > 0)
	#   i is the index of the qubit to be moved to the front (-n <= i < n)
	# Output:
	#   V which is overwritten
	#   On input, V = |Q0 Q1 ... Qi-1 Qi Qi+1 ... Qn-1>
	#   On output, V = |Qi Q0 Q1 ... Qi-1 Qi+1 ... Qn-1>

	n = qubitLength(V)

	assert -n <= i < n

	if i < 0:
		i += n

	n1 = n - 1
	ni = n1 - i
	hi_mask = (((1 << n1) - 1) >> ni) << ni
	lo_mask = (1 << ni) - 1

	newState = [V[((j & hi_mask) << 1) + (j & lo_mask) + ((j >> n1) << ni)] for j in xrange(1 << n)]

	for j, e in enumerate(newState):
		V[j] = e

	return V

def prepareBell(initialState=0):
	# Input:
	#   An optional integer between 0 and 3 representing the initial state
	#   of the two qubits to be entangled (default is 0)
	# Output:
	#   A state vector for two qubits entangled in one of the four Bell states
	#   (depending on the initial state):
	#
	#   Initial state                Bell state
	#   ------------------------     ----------
	#   00 => q1=[1,0], q2=[1,0] --> 00+11 => (|00> + |11>)/sqrt(2)
	#   01 => q1=[1,0], q2=[0,1] --> 01+10 => (|01> + |10>)/sqrt(2)
	#   10 => q1=[0,1], q2=[1,0] --> 00-11 => (|00> - |11>)/sqrt(2)
	#   11 => q1=[0,1], q2=[0,1] --> 01-10 => (|01> - |10>)/sqrt(2)

	assert 0 <= initialState <= 3

	q1 = [1, 0] if (initialState & 2) == 0 else [0, 1]
	q2 = [1, 0] if (initialState & 1) == 0 else [0, 1]

	changeState(HadamardGate, q1)
	return changeState(ControlledNotGate, combineStates(q1, q2))

def encodeBell(bit1, bit2, V):
	# Input:
	#   bit1 is a classical bit (0 or 1)
	#   bit2 is a classical bit (0 or 1)
	#   V is a state vector for one or more qubits
	# Output:
	#   None, but the input state vector V is overwritten

	assert bit1 == 0 or bit1 == 1
	assert bit2 == 0 or bit2 == 1

	vlen = validState(V)

	if bit1 == 0:
		if bit2 == 0:
			U = IdentityMatrix
		else:
			U = XGate
	else:
		if bit2 == 0:
			U = ZGate
		else:
			U = multiplyMatrixByMatrix(ZGate, XGate)

	while vlen > 2:
		U = combineTransforms(U, IdentityMatrix)
		vlen >>= 1

	changeState(U, V)

def measureQubit(V):
	# Input:
	#   V is a state vector for one or more qubits
	# Output:
	#   A classical bit (0 or 1) representing a measurement on the first qubit
	#   The input state vector V is overwritten

	vlen = validState(V)

	# The probability that the measurement will be 0 is the sum of the squares
	# of the absolute values of the first len(v)/2 elements of the state vector

	prob0 = sum([abs(e) * abs(e) for e in V[:vlen/2]])
	prob1 = sum([abs(e) * abs(e) for e in V[vlen/2:]])

	assert 1.0 - (prob0 + prob1) < 0.000000000001

	if random.random() < prob0:
		measurement = 0
		V[vlen/2:] = []
	else:
		measurement = 1
		V[:vlen/2] = []

	norm = math.sqrt(sum([abs(e) * abs(e) for e in V]))
	for i, e in enumerate(V):
		V[i] = e / norm

	return measurement

def measureBell(V):
	# Input:
	#   V is a state vector for two or more qubits
	# Output:
	#   The tuple (b1, b2) representing a measurement on the first two qubits
	#   The input state vector V is overwritten

	vlen = validState(V, 4)

	# Apply a controlled NOT gate on the first two qubits
	# followed by a Hadamard gate on the first qubit

	U = combineTransforms(HadamardGate, IdentityMatrix)
	U = multiplyMatrixByMatrix(U, ControlledNotGate)

	while vlen > 4:
		U = combineTransforms(U, IdentityMatrix)
		vlen >>= 1

	changeState(U, V)

	b1 = measureQubit(V) # Measure the first qubit
	b2 = measureQubit(V) # Measure the second qubit

	return b1, b2

def superdenseCoding(a1, a2):
	print '---- superdenseCoding() ----'
	print a1, a2
	qAB = prepareBell()
	encodeBell(a1, a2, qAB) # Alice
	(b1, b2) = measureBell(qAB) # Bob
	print b1, b2

def quantumTeleportation(qC):
	print '---- quantumTeleportation() ----'
	printQubit(qC)

	qAB = prepareBell()
	qABC = combineStates(qAB, qC)

	printQubit(qABC)
	reorderState(qABC, 2) # reorder AB12...n -> 1AB2...n
	printQubit(qABC)

	(b1, b2) = measureBell(qABC) # Alice
	encodeBell(b1, b2, qABC) # Bob

	print b1, b2
	printQubit(qABC)
	return qABC

def teleportEntangled():
	print '---- teleportEntangled() ----'

	q12 = prepareBell()
	qB2 = quantumTeleportation(q12)

	# Now qB is entangled with q2
	b1 = measureQubit(qB2) # Measure qB
	b2 = measureQubit(qB2) # Measure q2
	print b1, b2

def testRandomness():
	count = [0, 0]
	for i in xrange(0, 10000):
		v = prepareBell()
		b1 = measureQubit(v)
		b2 = measureQubit(v)
		assert b1 == b2
		count[b1] += 1
	print '0:', count[0]
	print '1:', count[1]

if __name__ == '__main__':
	superdenseCoding(1, 1)
	quantumTeleportation([0.6, 0.8])
	teleportEntangled()
