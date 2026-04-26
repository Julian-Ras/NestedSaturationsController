"""
This script is based on the article
Nested saturation with guaranteed real poles
https://ieeexplore.ieee.org/document/1239062
"""

import numpy as np
import matplotlib.pyplot as plt

# --- Saturation function ---
def sat(value, limit):
    return np.clip(value, -limit, limit)

# --- System matrices ---
A = np.array([[0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1],
              [0, 0, 0, 0]])

B = np.array([0, 0, 0, 9.8])

# --- Parameters from paper ---
p = 2 * np.array([0.5, 1, 2, 3])     # poles
b = 100 * np.array([0.625, 1.25, 2.5, 5])  # saturation limits

# --- Transformation matrix Tx ---
Tx = np.array([
    [p[0]*p[1]*p[2]*p[3],
     (p[0]*p[2] + p[1]*p[2] + p[0]*p[1]) * p[3],
     (p[0] + p[1] + p[2]) * p[3],
     p[3]],

    [0,
     p[0]*p[1]*p[2],
     (p[0] + p[1]) * p[2],
     p[2]],

    [0, 0,
     p[0]*p[1],
     p[1]],

    [0, 0, 0,
     p[0]]
])

# --- Simulation parameters ---
dt = 0.001
t_final = 50
s = np.arange(0, t_final + dt, dt)

# --- State variables ---
x = np.zeros((4, len(s)))
e = np.zeros((4, len(s)))
r = np.zeros((4, len(s)))
zx = np.zeros((4, len(s)))
u = np.zeros(len(s))

# Initial condition (already zero, but explicit)
x[:, 0] = np.array([0, 0, 0, 0])

# --- Simulation loop ---
for i in range(len(s) - 2):
    t = s[i]

    # Reference
    r[:, i] = np.array([
        4 + np.sin(0.5 * t),
        0.5 * np.cos(0.5 * t),
        -0.25 * np.sin(0.5 * t),
        -0.125 * np.cos(0.5 * t)
    ])

    # Error
    e[:, i] = x[:, i] - r[:, i]

    # 4th derivative reference
    xpppp = 0.0625 * np.sin(0.5 * t)

    # Transformation
    zx[:, i] = Tx @ e[:, i]

    # --- Nested saturation ---
    inner1 = sat(zx[0, i], b[0])
    inner2 = sat(zx[1, i] + inner1, b[1])
    inner3 = sat(zx[2, i] + inner2, b[2])
    inner4 = sat(zx[3, i] + inner3, b[3])

    u[i] = xpppp - inner4 / 9.8

    # System dynamics (Euler integration)
    x_dot = A @ x[:, i] + B * u[i] - np.array([0, 0, 0, xpppp])
    x[:, i + 1] = x[:, i] + dt * x_dot

# --- Plot 1: Position ---
plt.figure()
plt.plot(s, x[0, :], label='x1')
plt.plot(s, r[0, :], label='x1_ref')
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.legend()
plt.grid()

# --- Plot 2: Errors ---
plt.figure()
for i in range(4):
    plt.plot(s, e[i, :], label=f'e{i+1}')
plt.xlabel("Time")
plt.ylabel("Errors")
plt.legend()
plt.grid()

# --- Plot 3: Control input ---
plt.figure()
plt.plot(s, u, label='u')
plt.xlabel("Time")
plt.ylabel("Control input")
plt.grid()

plt.show()