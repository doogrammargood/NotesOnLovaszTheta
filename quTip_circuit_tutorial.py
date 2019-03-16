#import IPython
from IPython.display import Image
from numpy import pi
from qutip import *
#print cphase(pi/2)
qc3 = QubitCircuit(3)
qc3.add_gate("CNOT", 1, 0)
qc3.add_gate("RX", 0, None, pi/2, r"\pi/2")
qc3.add_gate("RY", 1, None, pi/2, r"\pi/2")
qc3.add_gate("RZ", 2, None, pi/2, r"\pi/2")
qc3.add_gate("ISWAP", [1, 2])
qc3.png
#print Image(filename='images/cphase.png')
# qc3 = QubitCircuit(3)
# qc3.add_gate("CNOT", 1, 0)
# qc3.add_gate("RX", 0, None, pi/2, r"\pi/2")
# qc3.add_gate("RY", 1, None, pi/2, r"\pi/2")
# qc3.add_gate("RZ", 2, None, pi/2, r"\pi/2")
# qc3.add_gate("ISWAP", [1, 2])
# IPython.core.display.DisplayObject(qc3.png)
