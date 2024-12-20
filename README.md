# Orr-Sommerfeld-Equation-Solver
Part of my Term Paper for the course **AM6515-Boundary Layer Stability** at **IIT MADRAS**.
This repository contains a **MATLAB** based code to solve the Orr-Sommerfeld Equation. 
The file *blasius_shooting.m* contains a method to solve the **Blasius profile** by a shooting method and RK-4 intetegration (*rk4.m*). It then solves the Orr-Sommerfeld equation by using a central difference scheme to estimate the second and fourth order derivatives which is used for approximating the Boundary conditions for the system of equations.
The file *plane_poiseulle.m* contains a code for solving the Orr-Sommerfeld equation for a **Plane-Posieulle flow profile** by Spectral Methods using the Chebyshev Polynomial expansions. 
Follow the references for gaining a deeper understanding of how the code works. They helped me a lot.

## The Orr-Sommerfeld Equation 
- The Orr-Sommerfeld equation is given by:
  
$$(\overline{U} - c)\left(\frac{d^2}{dy^2} - k^2\right)\hat{\psi} - \hat{\psi}\frac{d^2 \overline{U}}{dy^2} = \frac{1}{ikR_e}\left(\frac{d^2}{dy^2} - k^2\right)^2\hat{\psi}$$

It is derived from the governing equations i.e, the continuity and momentum equations after non-dimensionalising with relevant spatial and temporal scales. $$\Psi$$ is the stream function.
- A normal mode solution is adopted:

$$\psi = \frac{\hat{\psi(y)}}{2}e^{i(kx - \omega t)} + C.C$$

The Orr-Sommerfeld equation after adopting this is:

$$(\overline{U} - c)\left(\frac{d^2}{dy^2} - k^2\right)\hat{\psi} - \hat{\psi}\frac{d^2 \overline{U}}{dy^2} = \frac{1}{ikR_e}\left(\frac{d^2}{dy^2} - k^2\right)^2\hat{\psi}$$

## Numerical Techniques
### Shooting Method
## Orr-Sommerfeld Equation Solution

**Matrix Formulation:**  
   - Rewrite the O.S equation as:  
     
     $$\Phi' = A(y)\Phi$$
     
     with:  

  ```math
\Phi = \begin{pmatrix} \hat{\psi} & \frac{d\hat{\psi}}{dy} & \frac{d^2\hat{\psi}}{dy^2} & \frac{d^3\hat{\psi}}{dy^3} \end{pmatrix}, \quad A(y) = \begin{bmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ a & 0 & b & 0 \end{bmatrix}
```

  
     

2. **General Solution:**  
   - The general solution is:  
     
     $$\Phi(y) = \gamma_1\Phi^1(y) + \gamma_2\Phi^2(y) + \gamma_3\Phi^3(y) + \gamma_4\Phi^4(y)$$

3. **Initial Conditions:**  
   - At $$y = y_1$$  
   
     $$\Phi^1 = (0, 0, 0, 1)^T, \, \Phi^2 = (0, 0, 1, 0)^T, \, \Phi^3 = (0, 1, 0, 0)^T, \, \Phi^4 = (1, 0, 0, 0)^T$$

4. **Boundary Conditions:**  
   - Enforce physical boundary conditions:  
   
     $$\hat{\psi}(y_1) = 0, \quad \frac{d\hat{\psi}}{dy}(y_1) = 0$$ 
     This reduces the solution to:  
   
     $$\Phi(y) = \gamma_1\Phi^1(y) + \gamma_2\Phi^2(y)$$

5. **Numerical Integration:**  
   - Use RK-4 to solve for $$\Phi^1(y)$$  and $$\Phi^2(y)$$

6. **Eigenvalue Problem:**  
   - At $$y = y_2$$, impose:  
     
     $$\gamma_1\Phi^1_1(y_2) + \gamma_2\Phi^2_1(y_2) = 0, \quad \gamma_1\Phi^1_2(y_2) + \gamma_2\Phi^2_2(y_2) = 0$$

7. **Solve Determinant:**  
   - For non-trivial $$\gamma_1, \gamma_2$$ , solve:  

     $$L = \Phi^1_1(y_2)\Phi^2_2(y_2) - \Phi^1_2(y_2)\Phi^2_1(y_2) = 0$$

8. **Alternative Approach:**  
   - Use central differencing to approximate derivatives and simplify numerical implementation.



References:
1) [**Temporal and Spatial Stability Analysis of the Orr-Sommerfeld Equation**](https://www.cdsimpson.net/2015/04/temporal-and-spatial-stability-analysis.html#:~:text=This%20is%20a%20nonlinear%20ordinary,opposite%20boundary%20conditions%20are%20met).
2) [**Chebyshev collocation code for solving two phase Orr-Sommerfeld eigenvalue problem**](https://in.mathworks.com/matlabcentral/fileexchange/48862-chebyshev-collocation-code-for-solving-two-phase-orr-sommerfeld-eigenvalue-problem).
3) **Peter J. Schmid and Dan S. Henningson**, *Stability and Transition in Shear Flows*, *Springer*, 2001
