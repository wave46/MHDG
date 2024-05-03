$\kappa\text{-}\epsilon$

$$\begin{align}
    \frac{\partial\kappa}{\partial t} - \overbrace{\nabla_\perp \cdot (D\nabla_\perp \kappa)}^{\text{Diffusion}} + \overbrace{\nabla \cdot (\kappa u)}^{\text{Convection}}  &= \gamma \kappa - \frac{\kappa²}{D_\omega} - \epsilon \\
    \frac{\partial\epsilon}{\partial t} - \nabla_\perp \cdot (D\nabla_\perp \epsilon) + \nabla \cdot (\epsilon u) &= \gamma \epsilon - V \epsilon² \kappa^{-\frac{3}{2}}
\end{align}$$

where : 

$$\begin{align}
    D_\omega &= \frac{\kappa_{\max}}{\gamma} \\
    V &= c_s \rho_0 \frac{V_0}{R_0} \sqrt{\tau_\parallel \gamma} \\
    \gamma &= \dfrac{-A + \Re \left( \sqrt{A² - a_n a_\phi - \dfrac{b_n b_\phi}{C_\parallel / \Omega_i} }\right)}{\tau} \\
    A &= \frac{a_n + a_\phi }{2} \\
    a_n &= \frac{D_\perp}{D_B d_\perp}  + \sqrt{\frac{C_\parallel}{\Omega_i}}\\
    a_\phi &= \frac{\nu_\perp}{D_B d_\perp}  + d_\perp\\
    b_n &= \frac{\rho_0}{n_B} \left(2 d_\perp \right)^{-1/2} \Big( \frac{\partial n}{\partial R} -  \frac{\partial n}{\partial z}\Big)
    + i\left(\frac{C_\parallel}{\Omega_i}\right)^{3/4}\\
    b_\phi &= -\sqrt{2 d_\perp}\frac{\rho_0}{B} \frac{\partial B}{\partial R} + i\delta d_\perp \left(\frac{C_\parallel}{\Omega_i}\right)^{3/4}\\
\end{align}$$

## Dependencies
```
sudo apt-get install gfortran openmpi-bin openmpi-common openmpi-doc libopenmpi-dev libblas-dev liblapack-dev libxt-dev libhdf5-serial-dev
```
## Building
```zsh
cd lib
source Make.inc/init_vars_libs.sh
make
```
