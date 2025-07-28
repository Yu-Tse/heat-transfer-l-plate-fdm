# Finite‑Difference Analysis of an L‑Shaped Plate

<p align="center">
  <img src="docs/assets/problem_setup.png" width="450" alt="Problem setup: geometry, BCs & heat flux"/>
</p>

> **Author**
> • Code & report — Project members (see `CONTRIBUTORS.md`)
> • **Analytical derivation of the 17‑node finite‑difference formulas** — *Yu‑Chih Chi (紀詠喆)*, Dept. of Mechanical Engineering, **National Chung Cheng University**

---

## 1. Problem Statement

| Symbol     | Value         | Description                                         |
| ---------- | ------------- | --------------------------------------------------- |
| `L`        | 0.20 m        | Overall plate length                                |
| `T₀`       | 300 K         | Uniform initial temperature                         |
| `Tₛ`       | 400 K         | Prescribed temperature on the **west** side (green) |
| `T_∞`      | 300 K         | Ambient temperature for convection boundaries       |
| `h`        | 20 W m⁻² K⁻¹  | Convection coefficient                              |
| $\dot q''$ | 2 000 W m⁻²   | Uniform surface heat flux on the north edge         |
| `k`        | 15 W m⁻¹ K⁻¹  | Thermal conductivity                                |
| Δx = Δy    | 5 mm (≤ 5 mm) | Uniform grid size                                   |

Boundary conditions are summarised in the figure above:

* **Green** — constant temperature `Tₛ`
* **Red arrows** — imposed surface flux $\dot q''$
* **Blue clouds** — convection (`h, T_∞`)
* **Hatched faces** — adiabatic (`∂T/∂n = 0`)

The lower‑right quarter (0.1 m × 0.1 m) is removed, giving an L‑shaped domain that is discretised into **41 × 41 nodes** (1 681 in total; 400 lie inside the cut‑out and are therefore ignored).

---

## 2. Governing Equation

The transient 2‑D heat‑conduction equation with a volumetric heat source $\dot q'''$ reads

$$
\rho C_p\,\frac{\partial T}{\partial t}
\;=\;k\,\left(\frac{\partial^2 T}{\partial x^2}+\frac{\partial^2 T}{\partial y^2}\right)
\;\; +\; \dot q''',\qquad (1)
$$

where material properties are constant and isotropic.  In this project the source term is re‑defined as

$\dot q'''\;=\;\frac{2\,\dot q''}{\Delta x}\,\tag{2}$

to mimic the *effective* volumetric generation used in the MATLAB template.

### 2.1 Non‑dimensional form

Introducing the Fourier number $F_0 = \dfrac{k\,\Delta t}{\rho C_p\,\Delta x^2}$, Eq. (1) is cast into the explicit finite‑difference stencil

$$
T_{i,j}^{n+1}
= T_{i,j}^n
+ F_0\bigl(T_{i+1,j}^n+T_{i-1,j}^n+T_{i,j+1}^n+T_{i,j-1}^n-4T_{i,j}^n\bigr)
+ \dfrac{\Delta t}{\rho C_p}\,\dot q'''.\tag{3}
$$

The stability criterion for a square explicit mesh is $F_0 \le 0.25$; here $F_0 \approx 0.097<0.25$, hence the scheme is stable.

---

## 3. **17‑Node Formulation** (credit: *Yu‑Chih Chi*)

Because the L‑domain breaks regular connectivity, interior nodes fall into **17 distinct topological types**.  Chi’s derivation applies energy balance to a control volume around each type, leading to the following generic algebraic form

$a_P T_P = \sum a_{nb}\,T_{nb} + b\tag{4}$

where `P` is the current node, `nb` its neighbours (`W,E,N,S`), and coefficients `a` & source `b` depend on:

* Location (interior, boundary, cut‑out interface)
* Boundary condition (Dirichlet, Neumann/adiabatic, convection)
* Local heat flux / generation

The full derivation is provided in [`docs/derivation.pdf`](docs/derivation.pdf) ; a condensed table is reproduced below for quick reference.

| Node ID | Location / BC          | Finite‑difference coefficients                               |
| ------- | ---------------------- | ------------------------------------------------------------ |
| 1       | West + Dirichlet       | $T_P = T_s$                                                  |
| 2       | West interior          | `a_P = 1–4F₀`, `a_W = F₀`, `a_E = F₀`, `a_N = 2F₀`, `b = N₁` |
| …       | …                      | …                                                            |
| 17      | East convection corner | `a_P = 1–4F₀–2BiF₀`, `a_W = 2F₀`, `a_S = 2F₀`, `b = 2BiF₀T∞` |

*(see full table in the PDF)*

---

## 4. Code Usage

```bash
# clone & install deps
$ git clone https://github.com/your‑org/l‑plate‑fdm.git
$ cd l‑plate‑fdm
$ pip install -r requirements.txt  # numpy, matplotlib, pillow

# run the solver (≈ 1 min on a laptop)
$ python src/l_plate_explicit.py
```

Outputs:

* `steady_state.png`   – static 3‑D surface with colour bar
* `steady_rotation.gif` – 360‑degree rotating view *(embedded below)*
* `temp_evolution.gif` – 2‑D temperature evolution animation

<p align="center">
  <img src="output/steady_rotation.gif" width="420" alt="Rotating 3‑D steady‑state temperature"/>
</p>

---

## 5. Results & Discussion

| Property                           | Value                              |
| ---------------------------------- | ---------------------------------- |
| Steady‑state peak temperature      | **≈ 480 K** at the top‑left corner |
| Time to steady state (‖ΔT‖<10⁻⁴ K) | \~ 3 900 s                         |
| Stability check                    | `F0 = 0.0968 < 0.25` ✔︎            |

Key observations:

* The cut‑out corner (node‑9) exhibits the steepest gradient because heat can only leave through two directions.
* Increasing `k` shortens the transient period, whereas increasing `ρ` or `C_p` delays it (thermal inertia).
* The explicit scheme is memory‑light (one vector) but bound by the stability limit; a Crank–Nicolson or ADI variant can achieve larger `Δt`.

See the [presentation](docs/HeatTransfer_Final.pptx) for a deeper parametric study.

---

## 6. Repository Layout

```
├── docs/
│   ├── assets/          #   problem_setup.png, node_table.png, …
│   └── derivation.pdf   # – full 17‑node coefficient derivation (Chi, CCU)
├── output/              #   generated PNG / GIF after running the solver
├── src/
│   └── l_plate_explicit.py     #   main Python implementation
├── requirements.txt
└── README.md
```

---

## 7. License

This project is released under the **MIT License**.  Please cite appropriately if you build upon the derivation by *Yu‑Chih Chi (2024)*.
