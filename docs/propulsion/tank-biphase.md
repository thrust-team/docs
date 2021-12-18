---
title: Biphase Model for Saturated Tank
author: Alberto Scompairn
---
# Biphase Model for Saturated Tank
The tank holds a self pressurized oxidizer, which is in saturation condition between liquid and vapor.

## Oxidizer Extraction Model

Suppose a tank is filled up to a certain point with liquid oxidizer and its vapor is saturated at a certain ambient temperature.

### Parameters
The parameters of the tank before the extraction are:

- $V_T$ Internal volume of the tank
- $m_O$ Mass of oxidizer stored in the tank
- $T_0$ Initial extraction temperature

Other variables can thus be described, such as:

- $V_O$ Volume occupied by liquid oxidizer
- $f$ Ullage volume percentage

The following solver needs as initial conditions only:

- $T_0$ Starting saturation temperature
- $f$ Vapor volume percentage (ullage)

and the following data for the chosen oxidizer.

### Oxidizer Data
Data for known properties of the oxidizer as a function of temperature are needed:

- $p_{Sat}(T)$ Saturation pressure for a given temperature
- $\rho_l(T)$ Density of the saturated liquid for a given temperature
- $\rho_v(T)$ Density of the saturated vapor for a given temperature
- $s_l(T)$ Specific entropy of the saturated liquid for a given temperature
- $s_v(T)$ Specific entropy of the saturated vapor for a given temperature

### Numerical solution
#### Time constraints
The extraction time $t_n$ can be considered adimensional

$$ t_n = 1$$

therefore the $n$ computation steps are

$$n = \frac{t_n}{\Delta t}$$

where $\Delta t$ is a reasonably low number (ex. $\Delta t = 5\cdot10^{-4}$).

Although the initial conditions are useful, some of them can be discarded by considering the problem adimensionalized with respect to the mass of oxidizer stored inside the tank.

#### Initial conditions
$$ m_0 = 1 $$

Therefore the liquid mass is
$$ m_{l,0} = \frac{m_0}{1 + f \cdot \frac{\rho_{v,0}}{\rho_{l,0}}} $$
And vapor mass is
$$ m_{v,0} = m_0 - m_{l,0}$$

The oxidizer liquid volume is simply
$$ V_{l,0} = \frac{m_{l,0}}{\rho_{l,0}} $$
which can be used to get total volume
$$ V_0 = V_{l,0} \cdot (1+f) $$
and the vapor volume
$$ V_{v,0} = V_0 - V_{l,0} $$


The physical quantities are regulated by the following expressions:

$$
\begin{align}
v_l(T) &= \frac1{\rho_l(T)} \\
v_v(T) &= \frac1{\rho_l(T)} \\
x(v, T) &= \frac{v - v_l(T)}{v_v(T) - v_l(T)} \\
s(x, T) &= (1 - x) \cdot s_l(T) + x \cdot s_v(T) \\
v(x, T) &= (1 - x) \cdot v_l(T) + x \cdot v_v(T) \\
\end{align}
$$

and the following are the initial conditions of the loop

$$
\begin{align}
T_0 &\\
p_0 &= p_{Sat}(T_0)\\
v_0 &= \frac{V_0}{m_0} \\
x_0 &= x(v_0,T_0) \\
s_0 &= s(x_0, T_0) \\
S_0 &= m_0 \cdot s_0 \\
\dot{m}_0 &= \frac{m_0}{t_n}
\end{align}
$$

#### Mass extraction
The mass that gets extracted at every computation step from the tank is computed as

$$
dm_i = \dot{m}_i \cdot dt
$$

and the remaining mass inside the tank is thus:

$$
m_{i+1} = m_i - dm_i
$$

Whether the extraction is done of liquid oxidizer or vapor oxidizer, the entropy $S$ of the mass left in the tank is computed as:

$$
\begin{cases}
S_{i+1} = S_i - dm_i \cdot s_l(T_i),& \quad \text{liquid extraction}\\
S_{i+1} = S_i - dm_i \cdot s_v(T_i),& \quad \text{vapor extraction}
\end{cases}
$$

And thus the specific entropy $s$ and specific volume $v$ are

$$
\begin{align}
s_{i+1} = \frac{S_{i+1}}{m_{i+1}} \\
v_{i+1} = \frac{V_0}{m_{i+1}}
\end{align}
$$


#### Vapor quality and temperature evaluation
The next goal is to compute the variation of vapor quality $x$ and temperature $T$ given the variation in specific volume $v$ and specific enthalpy $s$ after th extraction.

The functions $s(x,T)$ and $v(x,T)$ expressed previously can be differentiated with respect to $x$ and $T$ as:

$$
\begin{align}
dv = \frac{\partial v}{\partial x}\, dx + \frac{\partial v}{\partial T}\,dT \\
ds = \frac{\partial s}{\partial x}\, dx + \frac{\partial s}{\partial T}\,dT
\end{align}
$$

which in matrix form becomes

$$
\begin{Bmatrix}
    dv \\
    ds
\end{Bmatrix}
= 
\begin{bmatrix}
    \frac{\partial v}{\partial x} &
    \frac{\partial v}{\partial T}\\
    \frac{\partial s}{\partial x}&
    \frac{\partial s}{\partial T}
\end{bmatrix}
\begin{Bmatrix}
    dx \\
    dT
\end{Bmatrix}
$$

and inverting algebraically

$$
\begin{Bmatrix}
    dx \\
    dT
\end{Bmatrix}
= 
\begin{bmatrix}
    \frac{\partial v}{\partial x} &
    \frac{\partial v}{\partial T}\\
    \frac{\partial s}{\partial x}&
    \frac{\partial s}{\partial T}
\end{bmatrix}^{-1}
\begin{Bmatrix}
    dv \\
    ds
\end{Bmatrix}
$$

The problem is now finding the quantities inside the Jacobian matrix. The partial derivatives with respect to vapor quality $x$ are immediate:

$$
\begin{align}
\frac{\partial v}{\partial x} = v_v(T) - v_l(T) \\
\frac{\partial s}{\partial x} = s_v(T) - s_l(T)
\end{align}
$$

but the derivatives with respect to temperature $T$ must be computed numerically, in this case with central finite difference:

$$
\begin{align}
\frac{\partial v}{\partial T}
&= \frac12\varepsilon [v(x,T+\varepsilon)-v(x,T-\varepsilon)] \\
\frac{\partial s}{\partial T}
&= \frac12\varepsilon [s(x,T+\varepsilon)-s(x,T-\varepsilon)]
\end{align}
$$

Starting the evaluation from the current temperature $T_i$, the target temperature $T_{i+1}$ will be reached with successive approssimation of a progressively changing temperature $T_j$ until the residue $r$ is under the desired tolerance:

$$
r = \sqrt{\frac{dT_j}{T_j}^2+\frac{dx_j}{x_j}^2}\leq \text{tol}
$$

The guesses for temperature $T_j$ and vapor quality $x_j$ are updated as:

$$
\begin{align}
T_{j+1} = T_j + dT_j \\
x_{j+1} = x_j + dx_j
\end{align}
$$
 until the final temperature $T_{i+1}$ is the $n$-th guess temperature $T_j\mid_{j=n}$.

 In summary:

$$
\begin{align}
T_{i+1} - T_{i} = \sum_{\substack{T_j=T_i\\x_j = x_i}}^{r<tol} dT_j(x_j,T_j) \\
x_{i+1} - x_{i} = \sum_{\substack{T_j=T_i\\x_j = x_i}}^{r<tol} dx_j(x_j,T_j)
\end{align}
$$

where, using the formula for the inverse of a $2 \times 2$ matrix:

$$
\begin{Bmatrix}
    dx_j \\
    dT_j
\end{Bmatrix}
=
\frac{1}{\frac{\partial v}{\partial x}\frac{\partial s}{\partial T} - \frac{\partial s}{\partial x}\frac{\partial v}{\partial T}}
\begin{bmatrix}
    \frac{\partial s}{\partial x} &
    -\frac{\partial v}{\partial T}\\
    -\frac{\partial s}{\partial x}&
    \frac{\partial v}{\partial T}
\end{bmatrix}
\begin{Bmatrix}
    dv_j \\
    ds_j
\end{Bmatrix}
$$

with all partial derivatives as a function of the guess vapor quality and temperature $f(x_j,T_j)$
and

$$
\begin{Bmatrix}
    dv_j \\
    ds_j
\end{Bmatrix}
=
\begin{Bmatrix}
    v_{i+1} - v(x_j,T_j) \\
    s_{i+1} - v(x_j,T_j)
\end{Bmatrix}
$$

#### Ready for next step
Once the quantities for $T_{i+1}$ and $x_{i+1}$ are estimated, the new value for pressure $p_i$ and mass flow $\dot{m}_i$ are:

$$
\begin{align}
p_{i+1} = p_{Sat}(T_{i+1}) \\
\dot{m}_{i+1} = \frac{\dot{m}_0}{p_0} \cdot p_{i+1}
\end{align}
$$

The last equation assumes that the mass flow ratio and pressure is constant along the entire extraction.

### Conclusion
The extraction is over once the vapor quality $x$ reaches the value of $1$, which is when the tank contains only vapor and the system is not biphasic anymore.

