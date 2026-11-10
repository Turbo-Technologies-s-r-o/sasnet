# sasnet — theory manual

Every model in sasnet follows from the conservation laws (mass, momentum,
energy) or from a classical law in the public domain. No handbook tables or
proprietary data are used in the code. Constants that cannot be derived are
named, sourced and overridable; project-specific calibration lives outside
the package.

Notation: $p$ total pressure, $T$ total temperature, $\rho$ density,
$u$ mean velocity, $A$ flow area, $d_h = 4A/P$ hydraulic diameter,
$\dot m$ mass flow, $Re = \rho u d_h/\mu$, $\zeta$ loss coefficient referred to
the channel head $\tfrac12\rho u^2$.

## 1. Network equations

Nodes are plenums (cavities) with negligible velocity, so node pressure and
temperature are total values. Branches connect two nodes.

**Mass** (Kirchhoff's first law) at every internal node $i$:
$$\sum_j \dot m_{ji} = 0 .$$

**Energy.** Adiabatic branches conserve total enthalpy; at a node the
inflows mix: $c_p T_i \sum \dot m_{in} = \sum \dot m_{in} c_p T_{in}$
(heat input along a branch adds $c_p\,\Delta T$).

**Solution.** Unknowns are the internal-node pressures. The residuals are
solved by Newton's method; the Jacobian is assembled branch by branch
(each branch depends only on its two node pressures), steps are limited and
back-tracked. Temperatures are updated in an outer fixed-point loop.

## 2. Branch momentum equation (Newton → Euler → Bernoulli)

Newton's second law for a fluid particle, written along a streamline
(Euler 1757), integrated between the ends of a branch with friction and
mixing losses lumped into $\zeta$:
$$p_1 - p_2 = \zeta\,\frac{\rho_m u^2}{2} - \Delta p_{c},$$
where $\Delta p_c$ is the centrifugal term (sec. 8) and $\rho_m$ the density
at the mean pressure. This is Bernoulli's equation with losses. Using the
mean density is exact to $O(M^2)$; the optional factor $(1 + M^2/4)$ on the
head is the next term of the compressible expansion of the dynamic pressure
($p_0 - p = \tfrac12\rho u^2(1 + M^2/4 + \dots)$). For given pressures the
equation is solved for $\dot m = \rho_m u A$ by iteration, since $\zeta$
depends on $Re$.

## 3. Contraction, vena contracta (Kirchhoff, Weisbach)

A jet leaving a sharp edge separates and contracts to $C_c A$. For a 2D slot
fed from an infinite reservoir, free-streamline theory (Kirchhoff 1869)
gives exactly
$$C_c = \frac{\pi}{\pi + 2} = 0.611 .$$
For an axisymmetric hole and a finite area ratio $F = A/A_{up}$ the
classical empirical law of Weisbach (1855) is used:
$$C_c = 0.63 + 0.37\,F^3 .$$

## 4. Sudden expansion (Borda–Carnot) with a real velocity profile

Momentum balance on a control volume from the exit of a channel (area $A_1$,
non-uniform profile) to the fully mixed section downstream (area $A_2$),
with the upstream pressure acting on the step face:
$$p_1 A_2 + \beta\rho u_1^2 A_1 = p_2 A_2 + \rho u_2^2 A_2 ,$$
and the energy flux $\alpha\,\tfrac12\rho u_1^3 A_1$. With $F_2 = A_1/A_2$
the total-pressure loss referred to $\tfrac12\rho u_1^2$ is
$$\zeta_{exp} = \alpha - 2\beta F_2 + F_2^2 ,$$
which reduces to Borda–Carnot $(1 - F_2)^2$ for a uniform profile
($\alpha = \beta = 1$) and to $\alpha$ for an exit into a plenum.

Profile coefficients $\alpha = \overline{(u/U)^3}$, $\beta = \overline{(u/U)^2}$:

* laminar, exact Poiseuille profiles: pipe $\alpha = 2$, $\beta = 4/3$;
  plates / narrow annulus $\alpha = 54/35$, $\beta = 6/5$;
* turbulent, power law $u/U_{max} = (y/R)^{1/n}$, integrated analytically:
  pipe
  $$\alpha = \frac{(n+1)^3(2n+1)^3}{4n^4(n+3)(2n+3)},\qquad
    \beta = \frac{(n+1)(2n+1)^2}{4n^2(n+2)},$$
  plates
  $$\alpha = \frac{(n+1)^3}{n^2(n+3)},\qquad \beta = \frac{(n+1)^2}{n(n+2)};$$
  the exponent follows from the friction law by Prandtl's rule
  $n = 1/\sqrt{\lambda}$ (sec. 6).

The developed-profile coefficients are applied for $l/d_h \ge 10$; shorter
channels use $\alpha = \beta = 1$.

## 5. Inlet, short channels, thin orifices

**Inlet with reattachment.** Flow entering a channel from a larger section
contracts to $C_c A$ and re-expands to $A$ (Borda–Carnot between the vena
contracta and the channel):
$\zeta = (1/C_c - 1)^2 \approx 0.35$ for $C_c = 0.63$. Measured sharp-edged
flush inlets give $0.5$ (Weisbach), the difference being mixing in the
separated region. sasnet uses
$$\zeta_{in} = K_{in}\,(1 - F_1),\qquad K_{in} = 0.5\ \text{(sharp edge)},$$
$F_1 = A/A_{up}$; a rounded edge reduces $K_{in}$ linearly to 0.04 at
$r/d_h = 0.15$. $K_{in}$ is the main calibration constant for holes.

**Thin orifice** ($l/d \to 0$, no reattachment). The jet leaves with
velocity $u/C_c$ and its kinetic energy above the downstream mean is lost:
$$\zeta_{orif} = \left(\frac{1}{C_c} - F_2\right)^2 ,$$
e.g. $2.52$ from a plenum to a plenum.

**Short channel.** Between the two limits the separated jet reattaches at
about $L_r = 2\,d_h$. The loss is blended quadratically:
$$\zeta = \zeta_{long} + (\zeta_{orif} - \zeta_{long})\,(1 - l/(L_r d_h))^2,
\quad l < L_r d_h ,$$
where $\zeta_{long} = \zeta_{in} + \lambda l/d_h + \zeta_{exp}$.

## 6. Friction (Hagen–Poiseuille, Prandtl–Nikuradse)

$\Delta p = \lambda\,(l/d_h)\,\tfrac12\rho u^2$ (Darcy–Weisbach).

* **Laminar**, exact solutions of the Navier–Stokes equations for developed
  flow: pipe $\lambda = 64/Re$ (Hagen 1839, Poiseuille 1840); parallel plates
  and narrow annulus $\lambda = 96/Re$; rectangle of aspect ratio $a$ by the
  series solution, fitted as
  $\lambda Re = 96(1 - 1.3553a + 1.9467a^2 - 1.7012a^3 + 0.9564a^4 - 0.2537a^5)$.
* **Turbulent, smooth**: Prandtl's universal law, obtained by integrating the
  logarithmic wall law over the cross-section with Nikuradse's constants,
  $$\frac{1}{\sqrt\lambda} = 2.0\,\lg\!\left(Re\sqrt\lambda\right) - 0.8 .$$
* **Rough**: Colebrook's combination of the smooth and fully rough log laws,
  $1/\sqrt\lambda = -2\lg\!\left(\varepsilon/(3.7d) + 2.51/(Re\sqrt\lambda)\right)$.
* **Transition** $2300 < Re < 4000$: linear blend.
* Non-circular turbulent ducts use $d_h$.

## 7. Bends (Weisbach)

$$\zeta_{bend} = \left[0.131 + 0.163\,(d/R)^{3.5}\right]\frac{\theta}{90^\circ},$$
$R$ the centreline radius (Weisbach 1855).

## 8. Heat addition and rotation

**Heating** at constant area, low Mach (Rayleigh line). Momentum:
$p_1 - p_2 = \rho_2 u_2^2 - \rho_1 u_1^2 = \rho_1 u_1^2 (T_2/T_1 - 1)$, since
$\rho u$ is constant and $\rho \propto 1/T$. Referred to the mean head,
$$\zeta_T = \frac{4\,(T_2 - T_1)}{T_1 + T_2}.$$

**Rotation.** In a frame rotating with the channel at $\omega$, the
centrifugal body force per unit volume is $\rho\omega^2 r$. Integrating
Euler's equation along a channel from $r_1$ to $r_2$ gives the pumping term
$$\Delta p_c = \tfrac12\rho\,\omega^2\,(r_2^2 - r_1^2),$$
a pressure gain for outward flow.

**Swirl entry.** Non-swirling fluid entering a channel moving at $\omega r$
must acquire the tangential velocity; the relative tangential kinetic energy
$\tfrac12\rho(\omega r)^2$ is (partly) lost. With $k = \omega r/u$:
$$\zeta_{rot} = c_{rot}\,k^2,\qquad 0 \le c_{rot} \le 1 .$$
$c_{rot} = 1$ is the upper bound; the actual value depends on the geometry
and must be set by the user (default 0, flagged in reports).

## 9. Labyrinth seals (Martin, jet carry-over, choking)

**Martin (1908).** $N$ identical throttles in series, isothermal ideal gas,
each throttle accelerating the gas from rest with all kinetic energy
dissipated in the following chamber. Summing $p_{k-1}^2 - p_k^2$ over the
throttles and using the continuity equation gives
$$\dot m = C_d\,\kappa\,A\,\frac{p_0}{\sqrt{RT_0}}
\sqrt{\frac{1 - (p_N/p_0)^2}{N - \ln(p_N/p_0)}} ,\qquad A = 2\pi r c .$$

**Carry-over.** In a straight-through seal part of the jet kinetic energy
passes over the next tooth. Assuming the jet spreads linearly with a small
half-angle (parameter 0.02), the fraction of the throttle cross-section
reached by the undisturbed jet is $X = (c/s)/(c/s + 0.02)$, and the effective
head increases by $1/(1 - X)$ for $N-1$ of the $N$ throttles:
$$\kappa = \left[1 - \frac{N-1}{N}\frac{c/s}{c/s + 0.02}\right]^{-1/2}$$
(Hodkinson 1939). Stepped seals: $\kappa = 1$.

**Choking.** Martin's formula assumes incompressible throttles and passes a
spurious maximum at low pressure ratio. sasnet treats the first $N-1$
throttles by Martin and the last one as an isentropic nozzle (sec. 10),
with the intermediate pressure found by bisection. This chokes the last
tooth physically and reduces to Martin at small pressure drops.

**Throttle coefficient.** $C_d$ of a sharp knife edge tends to the Kirchhoff
value 0.611; thick or rounded teeth give higher values. It is an input.

## 10. Isentropic nozzle (Saint-Venant–Wantzel)

From the energy equation $c_pT_0 = c_pT + u^2/2$ and the isentropic relation
$T/T_0 = (p/p_0)^{(\gamma-1)/\gamma}$:
$$\dot m = C_d A \frac{p_0}{\sqrt{RT_0}}
\sqrt{\frac{2\gamma}{\gamma-1}\left[\left(\frac{p}{p_0}\right)^{2/\gamma}
- \left(\frac{p}{p_0}\right)^{(\gamma+1)/\gamma}\right]} ,$$
maximal at $p/p_0 = (2/(\gamma+1))^{\gamma/(\gamma-1)} = 0.528$ for air
(Saint-Venant & Wantzel 1839); below it the flow is choked.

## 11. Rotor axial force

Each cavity face is an annulus $r..R$; the axial force is
$F = \pm\int_r^R p(\rho)\,2\pi\rho\,d\rho$. With uniform pressure
$F = \pm p\,\pi(R^2 - r^2)$. With a core rotating at $K\omega$, radial
equilibrium $dp/d\rho = \rho_g (K\omega)^2 \rho$ gives
$p(\rho) = p_{ref} + \tfrac12\rho_g(K\omega)^2(\rho^2 - r_{ref}^2)$, integrated
exactly. The rotor force is the sum over faces plus the blade force.

## Constants

| Constant | Value | Origin |
|---|---|---|
| $C_c$ (2D sharp slot) | $\pi/(\pi+2)$ | Kirchhoff 1869, exact |
| $C_c(F)$ | $0.63 + 0.37F^3$ | Weisbach 1855 |
| $K_{in}$ sharp inlet | 0.5 | Weisbach 1855 |
| reattachment length | $2\,d_h$ | separated-flow experiments, order of magnitude |
| laminar limit | $Re = 2300$ | Reynolds 1883 / Schiller 1922 |
| Prandtl law constants | 2.0, 0.8 | Prandtl 1935 / Nikuradse 1932 |
| jet spreading (labyrinth) | 0.02 | Hodkinson 1939 |
| critical pressure ratio | 0.528 ($\gamma = 1.4$) | exact |

## Validation policy

sasnet is checked against exact limits (tests) and against archived
project calculations and tests, which stay private. Agreement is judged as
trends relative to a known reference point. Calibrated constants
($K_{in}$, $C_d$ of seals, $c_{rot}$) belong to project data.
