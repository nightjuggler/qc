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

See the [`quantumFourierTransform()`](https://github.com/nightjuggler/qc/blob/92667c71095a66dbd80e3bbb51dd2cef9171b55b/qc.py#L450-L466) function.

## Acknowledgements

This library was inspired by Michael Nielsen's excellent video series [Quantum computing for the determined](http://michaelnielsen.org/blog/quantum-computing-for-the-determined/).

