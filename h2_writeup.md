# Homework 2: Darcy's Law, Permeability, and Hydraulic Conductivity

Viet M. Bui,
2026-09-06,
GLY 6826 - Hydrogeologic Modeling

Code source: [gh:vietbuiminh/hydro-modeling/porosity_and_rev.jl](https://github.com/vietbuiminh/hydro-modeling/blob/main/porosity_and_rev.jl)

> This matter because you can view my commits history
---

**1** A porous medium has an instrisic permeability of **1 darcy** [D]

a) The permeability of the medium in $cm^2$ 

$$
1[d] \approx 9.87 \cdot 10^{-13} [m^2] \approx 9.87 \cdot 10^{-9} [cm^2]
$$

b) Calculate the hydraulic conductivity of the medium for:

- water, with $v_w = 0.013 \frac{cm^2}{s}$

To calculate K we need to find the dynamic viscosity of water $\mu_w$, where density of water $\rho_w=1 [\frac{g}{cm^3}]$

$$
\mu_w = v_w \cdot \rho_w = 0.013 \left[\frac{g}{cm\cdot s}\right]
$$

$$
K_w = \frac{\rho_w g k}{\mu_w} = \frac{1 \cdot 981 \cdot 9.87\cdot 10^{-9}}{0.013} = 0.00074 \left[\frac{cm}{s}\right]
$$

$0.00074 \left[\frac{cm}{s}\right] \cdot 0.01 \left[\frac{m}{cm}\right]\cdot (60*60*24) \left[\frac{s}{day}\right]$
$= 0.63936 \left[\frac{m}{day}\right]$

- a viscous oil, $v_o = 1.8 \frac{cm^2}{s}$

find $\mu_o$ where $\rho_o = 0.956 [\frac{g}{cm^3}]$ referencing Bachaquero (Venezuela) heavy oil from [ZelenTech](https://www.zelentech.co/en/tools/oil-viscosity/) at 15°C 

$$
\mu_o = v_o \cdot \rho_o = 1.8 \cdot 0.956 = 1.7208 \left[\frac{g}{cm\cdot s}\right]
$$

$$
K_o = \frac{\rho_o g k}{\mu_o} = \frac{0.956 \cdot 981 \cdot 9.87\cdot 10^{-9}}{1.7208} \approx 0.00000538 \left[\frac{cm}{s}\right]
$$
$= 0.00464832 \left[\frac{m}{day}\right]$
