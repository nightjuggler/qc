#!/usr/bin/python
#
# qc.py - Quantum Computing Library
#
#         Various functions for mathematically simulating the quantum circuit model
#         of quantum computation including examples of superdense coding and quantum
#         teleportation.
#
#         February 13-18, 2016 - Pius Fischer
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

def changeLeadingState(U, V):
	vlen = validState(V)
	ulen = validMatrix(U)

	while vlen > ulen:
		U = combineTransforms(U, IdentityMatrix)
		ulen <<= 1

	changeState(U, V)

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

qubitStateMap = {}

class QubitState(object):
	def __init__(self, id, pa0, pa1):
		global qubitStateMap

		assert id not in qubitStateMap
		assert 1.0 - (abs(pa0)**2 + abs(pa1)**2) < 0.000000000001

		self.stateVector = [pa0, pa1]
		self.qubitNames = [id]

		qubitStateMap[id] = self

	def extend(self, otherState):
		if self is otherState:
			return

		self.stateVector = combineStates(self.stateVector, otherState.stateVector)
		self.qubitNames.extend(otherState.qubitNames)

		global qubitStateMap
		for id in otherState.qubitNames:
			qubitStateMap[id] = self

	def length(self):
		return len(self.qubitNames)

	def reorder(self, newOrder):
		if self.qubitNames[:len(newOrder)] == newOrder:
			return

		idMap = dict([(id, i) for i, id in enumerate(reversed(self.qubitNames))])

		newOrder = [(id, idMap.pop(id)) for id in newOrder]

		newOrder.extend([(id, idMap[id]) for id in self.qubitNames if id in idMap])

		bitShift = [((1 << bit), bit, shift) for bit, (id, shift) in enumerate(reversed(newOrder))]

		self.stateVector = [
			self.stateVector[
				sum([(((state & mask) >> bit) << shift) for mask, bit, shift in bitShift])
			]
			for state in xrange(len(self.stateVector))
		]

		self.qubitNames = [id for id, i in newOrder]

	def transform(self, unitaryMatrix):
		changeLeadingState(unitaryMatrix, self.stateVector)

	def measure(self, id):
		self.reorder([id])

		# The probability that the measurement will be 0 is the sum of the squares of
		# the absolute values of the first half of the elements of the state vector.

		V = self.stateVector
		vlen = len(V)

		prob0 = sum([abs(pa)**2 for pa in V[:vlen/2]])
		prob1 = sum([abs(pa)**2 for pa in V[vlen/2:]])

		assert 1.0 - (prob0 + prob1) < 0.000000000001

		if random.random() < prob0:
			measurement = 0
			V[vlen/2:] = []
			prob0, prob1 = 1, 0
		else:
			measurement = 1
			V[:vlen/2] = []
			prob0, prob1 = 0, 1

		norm = math.sqrt(sum([abs(pa)**2 for pa in V]))
		for i, pa in enumerate(V):
			V[i] = pa / norm

		del self.qubitNames[0]

		global qubitStateMap
		del qubitStateMap[id]
		QubitState(id, prob0, prob1)

		return measurement

	def printState(self, printed=None):
		n = len(self.qubitNames)

		print ','.join(self.qubitNames), '= ['
		for i, pa in enumerate(self.stateVector):
			print ('  {:0' + str(n) + 'b} -> {: }  p={}').format(i, pa, abs(pa)**2)
		print ']'

		if printed is not None:
			printed.extend(self.qubitNames)

def clearSystem():
	global qubitStateMap
	qubitStateMap.clear()

def createQubit(id, pa0, pa1):
	QubitState(id, pa0, pa1)

def removeQubit(id):
	global qubitStateMap
	assert qubitStateMap[id].length() == 1
	del qubitStateMap[id]

def printQubit(id):
	global qubitStateMap
	qubitStateMap[id].printState()

def printSystem():
	global qubitStateMap
	printed = []
	for id, state in qubitStateMap.iteritems():
		if id not in printed:
			state.printState(printed)

def applyGate(gate, *qubits):
	ulen = validMatrix(gate)
	qlen = 1 << len(qubits)
	assert ulen == qlen >= 2

	# Combine state vectors as necessary so all the qubits are in the same state vector:

	global qubitStateMap
	qState = qubitStateMap[qubits[0]]
	for id in qubits[1:]:
		qState.extend(qubitStateMap[id])

	qState.reorder(qubits)
	qState.transform(gate)

def measureQubit(id):
	global qubitStateMap
	return qubitStateMap[id].measure(id)

def prepareBell(q1, q2, initialState=0):
	# Input:
	#   q1 and q2 are the names to be given to the two entangled qubits.
	#
	#   initialState is an optional integer between 0 and 3 representing the
	#   initial (pre-entangled) state of the two qubits. The default is 0.
	#
	# Output: None
	#
	#   The two qubits q1 and q2 are created and entangled with each other in
	#   one of the four Bell states (depending on the initial state):
	#
	#   Initial state                Bell state
	#   ------------------------     ----------
	#   00 => q1=[1,0], q2=[1,0] --> 00+11 => (|00> + |11>)/sqrt(2)
	#   01 => q1=[1,0], q2=[0,1] --> 01+10 => (|01> + |10>)/sqrt(2)
	#   10 => q1=[0,1], q2=[1,0] --> 00-11 => (|00> - |11>)/sqrt(2)
	#   11 => q1=[0,1], q2=[0,1] --> 01-10 => (|01> - |10>)/sqrt(2)

	assert 0 <= initialState <= 3

	if (initialState & 2) == 0:
		createQubit(q1, 1, 0)
	else:
		createQubit(q1, 0, 1)

	if (initialState & 1) == 0:
		createQubit(q2, 1, 0)
	else:
		createQubit(q2, 0, 1)

	applyGate(HadamardGate, q1)
	applyGate(ControlledNotGate, q1, q2)

def encodeBell(bit1, bit2, qubit):
	# Input:
	#   "bit1" and "bit2" are classical bits (0 or 1) that determine which
	#   unitary matrix (quantum gate) to apply to the input "qubit".
	#
	# Output: None
	#
	assert bit1 == 0 or bit1 == 1
	assert bit2 == 0 or bit2 == 1

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

	applyGate(U, qubit)

def measureBell(q1, q2):
	# Input:
	#   q1 and q2 are the two qubits to be measured "in the Bell basis"
	# Output:
	#   The tuple (b1, b2):
	#   b1 is the result of the measurement on q1
	#   b2 is the result of the measurement on q2

	# Apply a controlled NOT gate on q1 and q2, with q1 as the control bit.
	# Then apply a Hadamard gate on q1.

	applyGate(ControlledNotGate, q1, q2)
	applyGate(HadamardGate, q1)

	b1 = measureQubit(q1)
	b2 = measureQubit(q2)

	return b1, b2

def sendSuperdense(a1, a2, senderQubit, receiverQubit):
	prepareBell(senderQubit, receiverQubit)
	encodeBell(a1, a2, senderQubit) # Alice
	return measureBell(senderQubit, receiverQubit) # Bob

def teleport(fromQubit, viaQubit, toQubit):
	prepareBell(viaQubit, toQubit)
	(b1, b2) = measureBell(fromQubit, viaQubit) # Alice
	encodeBell(b1, b2, toQubit) # Bob

def testSuperdenseCoding(a1, a2):
	#
	# Send two classical bits a1 and a2 via superdense coding using entangled qubits A and B
	#
	qA, qB = 'A', 'B'

	clearSystem()

	b1, b2 = sendSuperdense(a1, a2, qA, qB)

	assert a1 == b1 and a2 == b2

def testRandomness():
	qA, qB = 'A', 'B'
	count = [0, 0]
	for i in xrange(0, 10000):
		clearSystem()
		prepareBell(qA, qB)
		b1 = measureQubit(qA)
		b2 = measureQubit(qB)
		assert b1 == b2
		count[b1] += 1
	print '0:', count[0]
	print '1:', count[1]

if __name__ == '__main__':
	qA, qB, qC, qD = 'A', 'B', 'C', 'D'

	testSuperdenseCoding(0, 0)
	testSuperdenseCoding(0, 1)
	testSuperdenseCoding(1, 0)
	testSuperdenseCoding(1, 1)

	# Create a qubit C and teleport its state (via A) to qubit B
	clearSystem()
	createQubit(qC, 0.6, 0.8)
	teleport(qC, qA, qB)

	# Entangle A with B and teleport the state of A to D (via C)
	clearSystem()
	prepareBell(qA, qB)
	teleport(qA, qC, qD)

	# Now B is entangled with D
	b1 = measureQubit(qB)
	b2 = measureQubit(qD)
	assert b1 == b2
