# Finite‑Difference Analysis of an L‑Shaped Plate

> **Authors**
> • Code & report — Project members (see [`CONTRIBUTORS.md`](CONTRIBUTORS.md))
> • **Analytical derivation of the 17‑node finite‑difference formulas** — *Yu‑Chih Chi (紀詠喆)*, Dept. of Mechanical Engineering, **National Chung Cheng University**

---

## 1 · Problem statement

| Symbol                  | Value                              | Description                                         |
| ----------------------- | ---------------------------------- | --------------------------------------------------- |
| **L**                   | 0.20 m                             | Overall plate length                                |
| **T<sub>0</sub>**       | 300 K                              | Uniform initial temperature                         |
| **T<sub>s</sub>**       | 400 K                              | Prescribed temperature on the **west** side (green) |
| **T<sub>∞</sub>**       | 300 K                              | Ambient temperature for convection boundaries       |
| **h**                   | 20 W m<sup>−2</sup> K<sup>−1</sup> | Convection coefficient                              |
| \$\dot q''\$            | 2 000 W m<sup>−2</sup>             | Uniform heat flux on the north edge                 |
| **k**                   | 15 W m<sup>−1</sup> K<sup>−1</sup> | Thermal conductivity                                |
| \$\Delta x = \Delta y\$ | 5 mm                               | Uniform grid size                                   |

<img width="723" height="554" alt="image" src="https://github.com/user-attachments/assets/015c13ec-b1d3-467a-9df3-2ac06ce6f9a3" />

The lower‑right quarter (0.1 m × 0.1 m) is removed, giving an **L‑shaped domain** discretised into 41 × 41 nodes (1 681 total; 400 fall inside the cut‑out and are ignored).

---

## 2 · Governing equation

The transient 2‑D conduction equation with a volumetric source \$\dot q'''\$ is

$$
\rho C_p\,\frac{\partial T}{\partial t}
\;=\;k\,\Bigl(\tfrac{\partial^2 T}{\partial x^2}+\tfrac{\partial^2 T}{\partial y^2}\Bigr)
\; +\; \dot q'''. \tag{1}
$$

In this project we redefine the source as

$$
\dot q''' = \frac{2\,\dot q''}{\Delta x}, \tag{2}
$$

to match the MATLAB template.

### 2.1 Dimensionless explicit stencil

With the Fourier number \$F\_0 = \dfrac{k,\Delta t}{\rho C\_p,\Delta x^2}\$, Eq. (1) becomes

$$
T_{i,j}^{n+1} = T_{i,j}^n + F_0\Bigl(T_{i+1,j}^n + T_{i-1,j}^n + T_{i,j+1}^n + T_{i,j-1}^n - 4T_{i,j}^n\Bigr) + \frac{\Delta t}{\rho C_p}\,\dot q'''. \tag{3}
$$

The scheme is stable because \$F\_0 \approx 0.097 < 0.25\$.

---

## 3 · 17‑node formulation  *(Yu‑Chih Chi)*

Removing the quarter forces interior nodes into **17 connectivity types**.  Energy balance on each control volume yields the generic algebraic form

$$
a_P\,T_P = \sum_{nb} a_{nb}\,T_{nb} + b. \tag{4}
$$

Full derivation  → [`docs/derivation.pdf`](docs/derivation.pdf).

| Node ID | Location / BC          | Key coefficients                                                 |
| ------- | ---------------------- | ---------------------------------------------------------------- |
| 1       | West + Dirichlet       | \$T\_P = T\_s\$                                                  |
| 2       | West interior          | \$a\_P=1-4F\_0\$, \$a\_W=a\_E=F\_0\$, \$a\_N=2F\_0\$, \$b=N\_1\$ |
| …       | …                      | …                                                                |
| 17      | East convection corner | \$a\_P=1-4F\_0-2BiF\_0\$, \$a\_W=a\_S=2F\_0\$, \$b=2BiF\_0T\_∞\$ |

*(See PDF for the full table.)*

---

## 4 · Running the code

```bash
# clone & install deps
$ git clone https://github.com/your‑org/heat-transfer-l-plate-fdm.git
$ cd heat-transfer-l-plate-fdm
$ pip install -r requirements.txt   # numpy, matplotlib, pillow

# run the solver (≈ 1 min on a laptop)
$ python src/l_plate_explicit.py
```

Generated assets

| File                         | Preview                                     |
| ---------------------------- | ------------------------------------------- |
| `output/steady_state.png`    | 3‑D surface (static)                        |
| `output/steady_rotation.gif` | ![rotating gif](output/steady_rotation.gif) |
| `output/temp_evolution.gif`  | 2‑D transient animation                     |

---

## 5 · Results snapshot

| Metric                         | Value                      |
| ------------------------------ | -------------------------- |
| Peak steady temperature        | **≈ 480 K** (top‑left)     |
| Time to steady (‖ΔT‖ < 10⁻⁴ K) | \~3 900 s                  |
| Stability check                | \$F\_0 = 0.0968 < 0.25\$ ✅ |

---

## 6 · Repo layout

```
├── docs/
│   ├── assets/ …
│   └── derivation.pdf
├── output/              # generated PNG / GIF
├── src/
│   └── l_plate_explicit.py
├── CONTRIBUTORS.md
├── requirements.txt
└── README.md
```

---

## 7 · License

Released under the **MIT License**.  Cite *Yu‑Chih Chi (2024)* for the node‑derivation if reused.
