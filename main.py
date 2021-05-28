import tkinter as tk
from tkinter import Toplevel, messagebox
from tkinter.constants import E, LEFT, S, W

from threading import Thread

import math
from mpmath.libmp.libintmath import gcd
import sympy
import numpy as np

import pandas as pd
from matplotlib import pyplot as plt

from fractions import Fraction

from qiskit import BasicAer, QuantumCircuit, execute
from qiskit.visualization import plot_histogram, circuit_drawer


qasm_text = None

_a = None
_N = None

# Shor code copied from qiskit textbook

def is_prime(a):
    return a > 1 and all(a % i for i in np.arange(2, a))


def calculate():
    global _a
    global _N

    # check if not empty

    e_is_empty = False

    if len(e_textfield.get()) == 0:
        e_is_empty = True

    if e_is_empty:
        messagebox.showerror(title="Error", message="At least one of your entries is empty!")
        return

    # check all integer
    
    p = 3
    q = 5
    e = 0
    
    try:
        e = int(e_textfield.get())
    except Exception as e:
        messagebox.showerror(title="Error", message="At least one of your entries is not an integer!")
        return

    # check p is not q
    if p == q:
        messagebox.showerror(title="Error", message="p and q can't be the same number!")
        return

    # check p and q prime number
    if not is_prime(p):
        messagebox.showerror(title="Error", message="p is not a prime number!")
        return

    if not is_prime(q):
        messagebox.showerror(title="Error", message="q is not a prime number!")
        return

    # calculate N
    n = p * q

    # calculate φ(N)
    phi_n = (p - 1) * (q - 1)

    # check e "teilerfremd" to φ(N)
    gcd_e = math.gcd(phi_n, e)

    if gcd_e != 1:
        messagebox.showerror(title="Error", message="gcd(e, φ(N)) is not 1. φ(N) is " + str(phi_n) + "!")
        return

    # check 1 < e < φ(N)
    if e <= 1 or e >= phi_n:
        messagebox.showerror(title="Error", message="e hast to be in the following range 1 < e < φ(N)!")
        return

    # calculate d
    d = sympy.mod_inverse(e, phi_n)

    # print everything
    n_text_var.set("N = " + str(n))
    phi_text_var.set("φ(N) = " + str(phi_n))
    d_text_var.set("d = " + str(d))
    text_text_var.set("Private Key: (" + str(d) + ", " + str(n) + "); Public Key: (" + str(e) + ", " + str(n) + ")")

    a = None

    if len(a_text.get()):
        try:
            a = int(a_text.get())
        except Exception as e:
            messagebox.showerror(title="Error", message="'a' must be 2,7,8,11 or 13")
            return

    if a not in [2,7,8,11,13]:
        messagebox.showerror(title="Error", message="'a' must be 2,7,8,11 or 13")
        return

    _a = a
    _N = n

    Thread(target=crack_on_quantum_device, args=(n, a)).start()


def crack_on_quantum_device(n, a):
    calc_button_text.set("Cracking ...")
    
    backend = BasicAer.get_backend("qasm_simulator")

    n_count = 8

    # Create QuantumCircuit with n_count counting qubits
    # plus 4 qubits for U to act on
    qc = QuantumCircuit(n_count + 4, n_count)

    # Initialize counting qubits
    # in state |+>
    for q in range(n_count):
        qc.h(q)
        
    # And auxiliary register in state |1>
    qc.x(3+n_count)

    # Do controlled-U operations
    for q in range(n_count):
        qc.append(c_amod15(a, 2**q), 
                [q] + [i+n_count for i in range(4)])

    # Do inverse-QFT
    qc.append(qft_dagger(n_count), range(n_count))

    # Measure circuit
    qc.measure(range(n_count), range(n_count))


    job = execute(qc, backend=backend, shots=1024)
    result = job.result()

    calc_button_text.set("Calculate")

    rows, measured_phases = [], []
    for output in result.get_counts():
        decimal = int(output, 2)  # Convert (base 2) string to decimal
        phase = decimal/(2**n_count)  # Find corresponding eigenvalue
        measured_phases.append(phase)
        frac = Fraction(phase).limit_denominator(n)
        r = frac.denominator
        guess = [gcd(a**(r//2)-1, n), gcd(a**(r//2)+1, n)]

        # Add these values to the rows in our table:
        rows.append([f"{output}(bin) = {decimal:>3}(dec)", 
                    f"{decimal}/{2**n_count} = {phase:.2f}",
                    f"{frac.numerator}/{frac.denominator}",
                    r,
                    guess])
        
    # store the rows in a table
    headers=["Register Output", "Phase", "Fraction", "Guess for p", "Factor Guess"]
    df = pd.DataFrame(rows, columns=headers)
    result_string = df.__str__()

    open_new_result_window(result, result_string, n, qc)


def c_amod15(a, power):
    """Controlled multiplication by a mod 15"""
    if a not in [2,7,8,11,13]:
        raise ValueError("'a' must be 2,7,8,11 or 13")
    U = QuantumCircuit(4)        
    for iteration in range(power):
        if a in [2,13]:
            U.swap(0,1)
            U.swap(1,2)
            U.swap(2,3)
        if a in [7,8]:
            U.swap(2,3)
            U.swap(1,2)
            U.swap(0,1)
        if a == 11:
            U.swap(1,3)
            U.swap(0,2)
        if a in [7,11,13]:
            for q in range(4):
                U.x(q)
    U = U.to_gate()
    U.name = "%i^%i mod 15" % (a, power)
    c_U = U.control()
    return c_U


def qft_dagger(n):
    """n-qubit QFTdagger the first n qubits in circ"""
    qc = QuantumCircuit(n)
    # Don't forget the Swaps!
    for qubit in range(n//2):
        qc.swap(qubit, n-qubit-1)
    for j in range(n):
        for m in range(j):
            qc.cp(-np.pi/float(2**(j-m)), m, j)
        qc.h(j)
    qc.name = "QFT†"
    return qc


def open_new_result_window(device_result, result_string, n, circ):
    global global_result
    global global_circ

    global_result = device_result
    global_circ = circ

    result_window = Toplevel(root)
    result_window.title("Results")

    result_window_canvas = tk.Canvas(result_window, )
    result_window_canvas.grid(rowspan=2, columnspan=3)

    result_label = tk.Label(result_window, font=("Courier", 18), justify=LEFT, text=result_string)
    result_label.grid(row=0, column=0, columnspan=3)

    circ_button = tk.Button(result_window, text="Show Circuit", command=plot_circuit, font=("Arial", 18))
    circ_button.grid(row=1, column=0)

    period_button = tk.Button(result_window, text="Show Period", command=plot_period, font=("Arial", 18))
    period_button.grid(row=1, column=1)

    plot_button = tk.Button(result_window, text="Show Plot", command=plot_result, font=("Arial", 18))
    plot_button.grid(row=1, column=2)


def plot_period():
    global _a
    global _N

    x_bits = 4

    P = np.zeros((2**x_bits,1))
    for x in range(2**x_bits):
        P[x] = _a**x % _N
    
    fig, ax = plt.subplots()

    ax.plot(P)

    ax.set_xlabel('x')
    ax.set_ylabel(str(_a)+ '^x mod ' + str(_N))

    fig.show()

    fig.close()


def plot_result():
    global global_result

    fig = plot_histogram(global_result.get_counts())
    plt.show()

    plt.close()


def plot_circuit():
    global global_circ

    fig = circuit_drawer(global_circ, output='mpl', fold=-1)
    plt.show()

    plt.close()


def reset():
    e_textfield.delete(0, len(e_textfield.get()))
    a_text.delete(0, len(a_text.get()))

    n_text_var.set("N = p * q")
    phi_text_var.set("φ(N) = (p - 1) * (q - 1)")
    d_text_var.set("d: e * d ≡ 1 (mod φ(N))")
    text_text_var.set("Private Key: (d, N); Public Key: (e, N)")


root = tk.Tk()

root.title('RSA QuantaCrack')

canvas = tk.Canvas(root, )
canvas.grid(columnspan=6, rowspan=8)

p_label = tk.Label(root, text="p = 3", font=("Arial", 18))
p_label.grid(row=1, column=0, columnspan=2)

q_label = tk.Label(root, text="q = 5", font=("Arial", 18))
q_label.grid(row=1, column=2, columnspan=2)

e_label = tk.Label(root, text="e = ", font=("Arial", 18))
e_label.grid(row=1, column=4, sticky=E)

e_textfield = tk.Entry(root, width=2, font=("Arial", 18))
e_textfield.grid(row=1, column=5, sticky=W)

n_text_var = tk.StringVar()
n_label = tk.Label(root, textvariable=n_text_var, font=("Arial", 18))
n_label.grid(row=2, column=1, columnspan=4)

n_text_var.set("N = p * q")

phi_text_var = tk.StringVar()
phi_label = tk.Label(root, textvariable=phi_text_var, font=("Arial", 18))
phi_label.grid(row=3, column=1, columnspan=4)

phi_text_var.set("φ(N) = (p - 1) * (q - 1)")

d_text_var = tk.StringVar()
d_label = tk.Label(root, textvariable=d_text_var, font=("Arial", 18))
d_label.grid(row=4, column=1, columnspan=4)

d_text_var.set("d: e * d ≡ 1 (mod φ(N))")

text_text_var = tk.StringVar()
text_label = tk.Label(root, textvariable=text_text_var, font=("Arial", 18))
text_label.grid(row=5, column=0, columnspan=6)

text_text_var.set("Private Key: (d, N); Public Key: (e, N)")

a_label = tk.Label(root, text="a = ", font=("Arial", 18))
a_label.grid(row=6, column=2, sticky=E)

a_text = tk.Entry(root, width=2, font=("Arial", 18))
a_text.grid(row=6, column=3, sticky=W)

calc_button_text = tk.StringVar()
calc_button = tk.Button(root, textvariable=calc_button_text, command=calculate, font=("Arial", 18))
calc_button.grid(row=7, column=0, columnspan=2, sticky=S)

calc_button_text.set("Calculate")

reset_button = tk.Button(root, text="Reset", command=reset, font=("Arial", 18))
reset_button.grid(row=7, column=4, columnspan=2, sticky=S)

root.mainloop()
