#!/usr/local/bin/python3
#
# Various functions to test qc.py
#
import math
from qc import *

def testSuperdenseCoding(a1, a2):
	#
	# Send two classical bits a1 and a2 via superdense coding using entangled qubits A and B
	#
	qA, qB = 'A', 'B'
	clearSystem()

	b1, b2 = sendSuperdense(a1, a2, qA, qB)
	assert a1 == b1 and a2 == b2

def testTeleport1():
	#
	# Create a qubit C and teleport its state (via A) to qubit B
	#
	qA, qB, qC = 'A', 'B', 'C'
	clearSystem()

	createQubit(qC, 0.6, 0.8)
	printQubit(qC)

	teleportQubit(qC, qA, qB)
	printQubit(qB)

def testTeleport2():
	#
	# Entangle A with B and teleport the state of A to D (via C)
	#
	qA, qB, qC, qD = 'A', 'B', 'C', 'D'
	clearSystem()

	prepareBell(qA, qB)
	printQubit(qA)

	teleportQubit(qA, qC, qD)

	# Now D is entangled with B
	printQubit(qD)
	b1 = measureQubit(qD)
	printQubit(qD)
	printQubit(qB)
	b2 = measureQubit(qB)
	assert b1 == b2

def measure5050():
	qA, qB = 'A', 'B'
	prepareBell(qA, qB)
	b1 = measureQubit(qA)
	b2 = measureQubit(qB)
	assert b1 == b2
	return b1

def measure3664():
	createQubit('A', 0.6, 0.8)
	return measureQubit('A')

def measure23_13_23_41():
	qA, qB = 'A', 'B'
	createQubit(qA, 0.6, 0.8)
	createQubit(qB, 0.8, 0.6)
	applyGate(ControlledNotGate, qA, qB)
	b1 = measureQubit(qA)
	b2 = measureQubit(qB)
	return (b1 << 1) + b2

def testRandomness(getMeasurement):
	count = [0] * 16
	for i in range(10000):
		clearSystem()
		m = getMeasurement()
		count[m] += 1
	for m, c in enumerate(count):
		if c > 0:
			print('{:04b}: {}'.format(m, c))

def testQFT():
	clearSystem()

	x = qubitArray('x', 5)
	createQubit(x[0], 0.6, 0.8)
	createQubit(x[1], 0.8, 0.6)
	createQubit(x[2], 1.0, 0.0)
	prepareBell(x[3], x[4])
	applyGate(FourierTransform(1 << len(x)), *x)

	y = qubitArray('y', 5)
	createQubit(y[0], 0.6, 0.8)
	createQubit(y[1], 0.8, 0.6)
	createQubit(y[2], 1.0, 0.0)
	prepareBell(y[3], y[4])
	quantumFourierTransform(y)

	compareStateVectors(x[0], y[0])

def multiplyWith(F, gate, numBefore, numAfter):
	for i in range(numBefore):
		gate = combineTransforms(IdentityMatrix, gate)
	for i in range(numAfter):
		gate = combineTransforms(gate, IdentityMatrix)

	return multiplyMatrixByMatrix(gate, F)

def constructQFT(N):
	#
	# This function iteratively builds up the matrix F for the quantum Fourier transform
	# for vectors of length 2^N (N qubits), starting from the 1x1 identity matrix (which
	# is just F for N=0). F is multiplied by the matrices for the swap gate, the Hadamard
	# gate, and the controlled phase shift gates. Before multiplication with F, each of
	# these gate matrices, which act on either one qubit (the Hadamard) or two adjacent
	# qubits, is "padded" with 2x2 identity matrices for the qubits that are not acted on
	# by the gate. This padding is done by calling combineTransforms() (which returns the
	# Kronecker product of two matrices) in the above helper function.
	#
	# The matrix generated by this method should be equivalent to the one generated by the
	# FourierTransform function in qc.py (with argument 1<<N instead of N).
	#
	assert N >= 0
	#
	# At the beginning of each iteration i (0 <= i < N) of the outer loop, F is already
	# equal to the QFT matrix for i qubits. In other words, F == FourierTransform(1 << i)
	# (within some tolerance due to floating point rounding errors). Each iteration then
	# incorporates into F the operations needed to build the QFT matrix for one additional
	# qubit. Let's call this additional qubit x[i]. First, F is expanded by setting it to
	# the Kronecker product of itself and the 2x2 identity matrix. This doubles the number
	# of rows and columns in F, thus quadrupling the number of elements in F. Operating on
	# an additional qubit means operating on double the number of states. Each iteration j
	# (i >= j > 0) of the inner loop then incorporates one swap and one controlled phase
	# shift operation into F. This causes x[i] to be swapped with each previous input qubit
	# so that after i iterations of the inner loop, x[i] becomes y[0] (the i'th input qubit
	# becomes the 0'th output qubit). Also, after each swap, the controlled phase gate
	# R(pi/2^j) is applied to x[i] and the qubit with which it was just swapped (x[i-j]),
	# with x[i] as the control bit. Finally, a Hadamard transform is incorporated into F,
	# corresponding to the Hadamard gate being applied to x[i].

	F = [[1]]
	for i in range(N):
		F = combineTransforms(F, IdentityMatrix)

		for j in range(i, 0, -1):
			gate = ControlledGate(PhaseShiftGate(math.pi / (1 << j)))
			gate = multiplyMatrixByMatrix(gate, SwapGate)
			F = multiplyWith(F, gate, j - 1, i - j)

		F = multiplyWith(F, HadamardGate, 0, i)

	compareMatrices(F, FourierTransform(1 << N))

if __name__ == '__main__':
	testSuperdenseCoding(0, 0)
	testSuperdenseCoding(0, 1)
	testSuperdenseCoding(1, 0)
	testSuperdenseCoding(1, 1)

#	testTeleport1()
#	testTeleport2()

#	testRandomness(measure5050)
#	testRandomness(measure3664)
#	testRandomness(measure23_13_23_41)

	testQFT()
	constructQFT(5)
