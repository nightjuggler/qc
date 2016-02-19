#!/usr/bin/python
#
# Various functions to test qc.py
#
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
	b1 = measureQubit('A')
	b2 = measureQubit('B')
	return (b1 << 1) + b2

def testRandomness(getMeasurement):
	count = [0] * 16
	for i in xrange(0, 10000):
		clearSystem()
		m = getMeasurement()
		count[m] += 1
	for m, c in enumerate(count):
		if c > 0:
			print '{:04b}: {}'.format(m, c)

if __name__ == '__main__':
	testSuperdenseCoding(0, 0)
	testSuperdenseCoding(0, 1)
	testSuperdenseCoding(1, 0)
	testSuperdenseCoding(1, 1)

	testTeleport1()
	testTeleport2()

#	testRandomness(measure5050)
#	testRandomness(measure3664)
#	testRandomness(measure23_13_23_41)
