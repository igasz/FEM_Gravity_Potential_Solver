# Iga Szaflik 4.4

import numpy as np
import matplotlib.pyplot as plt
import sys

def solve():
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
    else:
        n = 10

    if n < 2:
        print("n musi być >= 2")
        return

    D = 3 #dziedzina
    G = 1
    h = D/n # b - a
    nodes = np.linspace(0, D, n+1) # pkt pomiędzy przedziałami

    # warunki dirichleta
    PHI_LEFT = 5
    PHI_RIGHT = 2

    # cx + d
    c = (PHI_RIGHT - PHI_LEFT) / D
    d = PHI_LEFT

    # funkcja gęstości
    def rho(x):
        if 0 <= x <= 1:
            return -10
        elif 1 < x <= 2:
            return 1
        elif 2 < x <= 3:
            return -10
        return 0

    gauss_pts = [-1 / np.sqrt(3), 1 / np.sqrt(3)]
    gauss_w = [1, 1]

    # całki liczone numerycznie
    def d_phi_thilde(x):
        return c

    def integral_L(node_type, el_idx):
        res = 0.0
        a = nodes[el_idx]

        for p, weight in zip(gauss_pts, gauss_w):
            x_mapped = a + (p + 1) * h / 2
            v = (1 - p) / 2 if node_type == "left" else (1 + p) / 2
            res += (4 * np.pi * G * rho(x_mapped) * v) * (h / 2) * weight
        return res

    def integral_L1(node_type, el_idx):
        res = 0.0
        a = nodes[el_idx]
        dv = -1 / h if node_type == "left" else 1 / h

        for p, weight in zip(gauss_pts, gauss_w):
            x_mapped = a + (p + 1) * h / 2
            res += (d_phi_thilde(x_mapped) * dv) * (h / 2) * weight
        return res

    # budowanie macierzy A*w = b
    A = np.zeros((n-1, n-1)) # macierz
    b = np.zeros(n-1) # wektor b

    for i in range(1, n):
        idx = i - 1

        A[idx, idx] = -2 / h
        if idx > 0: A[idx, idx - 1] = 1 / h
        if idx < (n-1) - 1: A[idx, idx + 1] = 1 / h

        b[idx] += integral_L("right", i - 1)
        b[idx] += integral_L("left", i)
        b[idx] += integral_L1("right", i - 1)
        b[idx] += integral_L1("left", i)

    # rozwiązanie układu równań
    w = np.linalg.solve(A, b)
    phi = np.zeros(n + 1)

    for i in range(n + 1):
        x = nodes[i]
        shift = c * x + d
        if i == 0 or i == n:
            phi[i] = shift
        else:
            phi[i] = w[i - 1] + shift

    # wypisywanie wyników obliczeń
    print(f"Liczba elementów (n): {n}")
    print(f"Szerokość elementu (h): {h:.4f}")
    print(f"{'Węzeł i':<10} | {'x_i':<10} | {'Phi(x_i)':<10}")
    print("-" * 33)
    for i in range(n + 1):
        print(f"{i:<10} | {nodes[i]:<10.4f} | {phi[i]:<10.4f}")

    # rysowanie wykresów
    plt.figure(figsize=(10, 6))
    plt.plot(nodes, phi, 'r-o', label=f'n={n}')
    plt.xlabel("Pozycja x")
    plt.ylabel("Wartość Phi(x)")
    plt.grid(True)
    plt.legend()
    plt.show()

solve()