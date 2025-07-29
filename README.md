# Finite‑Difference Analysis of an L‑Shaped Plate

> **Authors & roles**  
> • *Yu‑Chih Chi (紀詠喆)* — derived the 17‑node finite‑difference formulas and built the initial MATLAB prototype.  
> • Current maintainer — ported the solver to Python and added 3‑D/animated visualisation.  
> • This project was originally completed as a final report for the undergraduate Heat Transfer course at National Chung Cheng University (2023), by Yu-Chih Chi and [Yu-Tse Wu]. It is published here for reference and demonstration purposes.



## 1 · Problem statement

| Symbol       | Value        | Description                              |
| ------------ | ------------ | ---------------------------------------- |
| **L**        | 0.20 m       | Overall plate length                     |
| **T₀**       | 300 K        | Uniform initial temperature              |
| **Tₛ**       | 400 K        | Dirichlet boundary on west face          |
| **T∞**       | 300 K        | Ambient temperature for convection faces |
| **h**        | 20 W m⁻² K⁻¹ | Convection coefficient                   |
| \$\dot q''\$ | 2 000 W m⁻²  | Surface heat flux on the north edge      |
| **k**        | 15 W m⁻¹ K⁻¹ | Thermal conductivity                     |
| Δx = Δy      | 5 mm         | Uniform grid spacing                     |

<img width="723" height="554" alt="image" src="https://github.com/user-attachments/assets/4a631b46-1ef8-4bea-bf6a-756a91f41039" />


The lower‑right quarter (0.1 m × 0.1 m) is removed, leaving **1281 active nodes** on a 41 × 41 mesh.

---

## 2 · Governing equation (plain‑text blocks for safe rendering)

```text

(1)  ρ C_p ∂T/∂t  =  k ( ∂²T/∂x² + ∂²T/∂y² ) + q‴

      with  q‴  =  2 q″ / Δx          (2)

```

### 2.1 Explicit stencil (dimensionless form)

```text

Fourier number :  F₀ = k Δt / ( ρ C_p Δx² )

Update rule :
(3)  T_{i,j}^{n+1} = T_{i,j}^n + F₀ ( T_E + T_W + T_N + T_S – 4 T_P ) + ( Δt / ρ C_p ) q‴
                    
Stability :  F₀ ≈ 0.097  < 0.25   →  stable.

```

---

## 3 · 17‑node formulation *(Yu‑Chih Chi)*

```text

Generic energy balance for any node type :
(4)  a_P T_P  =  Σ a_nb T_nb  +  b

Seventeen distinct connectivity cases arise due to the L‑shaped cut‑out; their
coefficients (a_P, a_nb, b) follow Chi’s original derivation.

```

---

## Usage & results

Run `src/l_plate_explicit.py` (Python 3.8+) to generate:

* `output/steady_state.png` — static 3‑D temperature surface
* `output/steady_rotation.gif` — 360° rotating view
* `output/temp_evolution.gif` — 2‑D transient animation

Sample output:

![temp_evolution](https://github.com/user-attachments/assets/4d284041-495e-4094-977d-0965f42ca3ba)

---

Released under the **MIT License**.  If you reuse the 17‑node derivation, please cite *Yu‑Chih Chi, National Chung Cheng University (2023)*.

## 🙋‍♂️ Author

**Yu-Tse Wu** (吳雨澤)  
*Master’s Student at the Institute of Innovation and Semiconductor Manufacturing, National Sun Yat-sen University*

GitHub: [@Yu-Tse](https://github.com/Yu-Tse)
