# qc.py - Quantum Computing Library for Python

**qc.py** is a Python library for mathematically simulating the quantum circuit model of computation.

## Example 1 - [Superdense Coding](https://en.wikipedia.org/wiki/Superdense_coding)

Send two classical bits a1 and a2 from Alice to Bob via an entangled pair of qubits:

```python
from qc import *

qA, qB = 'A', 'B'             # Define the names to be used for the qubits.
a1, a2 = 1, 0                 # Initialize the classical bits to be sent.

prepareBell(qA, qB)           # Create the entangled pair qA and qB.
encodeBell(a1, a2, qA)        # Alice encodes a1 and a2 onto qA (which also affects qB).
b1, b2 = measureBell(qA, qB)  # Bob recovers b1 and b2 by Bell measurement of qA and qB.
```

This is implemented by the [`sendSuperdense()`](https://github.com/nightjuggler/qc/blob/92667c71095a66dbd80e3bbb51dd2cef9171b55b/qc.py#L440-L443) function.

## Example 2 - [Quantum Teleportation](https://en.wikipedia.org/wiki/Quantum_teleportation)

Transfer the state of a qubit qC from Alice to Bob via an entangled pair of qubits qA and qB and two classical bits b1 and b2:

```python
from qc import *

qA, qB, qC = 'A', 'B', 'C'    # Define the names to be used for the qubits.
createQubit(qC, 0.8, 0.6)     # Initialize the qubit to be teleported.

prepareBell(qA, qB)           # Create the entangled pair qA and qB.
b1, b2 = measureBell(qC, qA)  # Alice gets b1 and b2 by Bell measurement of qC and qA.
encodeBell(b1, b2, qB)        # Bob encodes b1 and b2 onto qB. Now qB is in the same
                              # state that qC was in before Alice's Bell measurement.
```

This is implemented by the [`teleportQubit()`](https://github.com/nightjuggler/qc/blob/92667c71095a66dbd80e3bbb51dd2cef9171b55b/qc.py#L445-L448) function.

## Example 3 - [Quantum Fourier Transform](https://en.wikipedia.org/wiki/Quantum_Fourier_transform)

The [`quantumFourierTransform()`](https://github.com/nightjuggler/qc/blob/92667c71095a66dbd80e3bbb51dd2cef9171b55b/qc.py#L450-L466) function implements the quantum Fourier transform on a list of qubits by applying Hadamard gates and controlled phase shift gates.

The [`FourierTransform(N)`](https://github.com/nightjuggler/qc/blob/5054083b953263a6613bca1267b11eb14e432e02/qc.py#L155-L160) function generates the NxN gate matrix for the quantum Fourier transform by computing powers of the primitive N'th [root of unity](https://en.wikipedia.org/wiki/Root_of_unity).

`quantumFourierTransform(x)` is equivalent to `applyGate(FourierTransform(1 << len(x)), *x)`

See [`testQFT()`](https://github.com/nightjuggler/qc/blob/5054083b953263a6613bca1267b11eb14e432e02/qc-test.py#L82-L99) in [qc-test.py](qc-test.py) for an example.

## Example 4 - printQubit() and measureQubit()

The state of a qubit can be displayed with `printQubit()`.
A measurement on a qubit can be simulated with `measureQubit()` which returns 0 or 1 with the probability given by the qubit's state vector.
The states of all qubits in the system can be displayed with `printSystem()`.

The following example is also implemented by the [`testTeleport2()`](https://github.com/nightjuggler/qc/blob/5054083b953263a6613bca1267b11eb14e432e02/qc-test.py#L31-L49) function in [qc-test.py](qc-test.py).

First, two qubits **A** and **B** are created and entangled with each other in the Bell state.

```
>>> from qc import *
>>> qA, qB, qC, qD = 'A', 'B', 'C', 'D'
>>> prepareBell(qA, qB)
>>> printQubit(qA)
A,B = [
  00 ->  0.707106781187  p=0.5
  01 ->  0  p=0
  10 ->  0  p=0
  11 ->  0.707106781187  p=0.5
]
```

The state of **A** is then teleported to **D** (via **C**). **D** is now entangled with **B**.

```
>>> teleportQubit(qA, qC, qD)
>>> printQubit(qD)
D,B = [
  00 ->  0.707106781187  p=0.5
  01 ->  0  p=0
  10 ->  0  p=0
  11 ->  0.707106781187  p=0.5
]
```

Since **A** and **C** were measured as part of the teleportation protocol, they are no longer in a superposition of states. In this instance, 0 was measured for **A**, and 1 was measured for **C**.

```
>>> printSystem()
A = [
  0 ->  1.0  p=1.0
  1 ->  0  p=0
]
C = [
  0 ->  0  p=0
  1 ->  1.0  p=1.0
]
D,B = [
  00 ->  0.707106781187  p=0.5
  01 ->  0  p=0
  10 ->  0  p=0
  11 ->  0.707106781187  p=0.5
]
```

If we now measure **D**, we have a 50/50 chance of getting 0 or 1. In this instance, we measured 0, which means that **B** must then also be 0.

```
>>> measureQubit(qD)
0
>>> measureQubit(qB)
0
>>> printSystem()
A = [
  0 ->  1.0  p=1.0
  1 ->  0  p=0
]
C = [
  0 ->  0  p=0
  1 ->  1.0  p=1.0
]
B = [
  0 ->  1.0  p=1.0
  1 ->  0  p=0
]
D = [
  0 ->  1.0  p=1.0
  1 ->  0  p=0
]
>>>
```

## Acknowledgements

This project was inspired by Michael Nielsen's excellent video series "[Quantum computing for the determined](http://michaelnielsen.org/blog/quantum-computing-for-the-determined/)".
