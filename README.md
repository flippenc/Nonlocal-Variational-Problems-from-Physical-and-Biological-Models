**Research Funding and Support Information:**

-   This research was funded by [NSF Grant DMS
    2306962](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2306962&HistoricalAwards=false)

-   Some of the numerical results for this research were obtained using
    the Teal Cluster from VCU’s High Performance Computing Facility.
    Information about the facility [can be found
    here](https://research.vcu.edu/cores/hprc/facilities/)

**Theory Explanation:**

The Hamiltonian of a particle interaction system with *N* particles is
given by
$$E(\boldsymbol{x_1}, \ldots, \boldsymbol{x_N}) = \sum\_{\substack{i,j = 1 \\ i \neq j}}^N K(\boldsymbol{x_i}-\boldsymbol{x_j})$$
where **x**<sub>**1**</sub>, …, **x**<sub>**N**</sub> are
two-dimensional vectors and *K* : ℝ → ℝ is the interaction kernel. To
determine the ground state of this system, we solve the following ODE
system
$$\frac{d\boldsymbol{x_i}}{dt} = -\frac{1}{N}\sum\_{\substack{i,j = 1 \\ i \neq j}}^N \nabla K(\boldsymbol{x_i}-\boldsymbol{x_j}).$$
The two types of initial conditions we considered were:

-   particles randomly placed in a ball in ℝ<sup>2</sup> with radius 0.5

    -   the files `parameterizedGenModel.m`,
        `parameterizedGenModelElliptical.m`,  
        `parameterizedGenModelGamma.m`,`parameterizedGenModelLambda.m`,  
        and `parameterizedGenModelNorm.m` use the ball initial
        conditions

-   particles randomly placed in two balls in ℝ<sup>2</sup> of radius
    *ε*/4 sitting next to each other (balls centered at
    ($\pm\sqrt{\varepsilon/4},0$))

    -   the file `twoCirclesParameterizedGenModel.m` uses the two balls
        initial condition

The types of kernels we considered are as follows.

-   the “spherical kernel" which uses parameters *α*, *β*, *γ*, *δ*,
    *ε*, *λ*

    -   the file `parameterizedSystemGrad.m` computes the gradient of
        this kernel

-   the “elliptical kernel" which uses parameters *α*, *β*, *γ*, *δ*,
    *ε*, *λ*, *a*, *b*

    -   the file `ellipticalParameterizedSystemGrad.m` computes the
        gradient of this kernel

-   the “elliptical *c*-norm kernel" which uses parameters *α*, *β*,
    *γ*, *δ*, *ε*, *λ*, *a*, *b*, *c*

    -   the file `normParameterizedSystemGrad.m` computes the gradient
        of this kernel

**Spherical Kernel**

Let *α*, *β*, *γ*, *δ*, *ε*, and *λ* be real numbers with 2 \> *α* ≥ 0,
*β* ≥ 0, *γ* ≥ 0, *δ* ≥ 1, 1 ≫ *ε* \> 0, and *λ* \> 0. The “spherical
kernel" is defined as
$$K(x,y) = -\frac{\delta}{4\sqrt{\pi}\varepsilon^{3/2}}\exp\left(-\frac{x^2+y^2}{4\varepsilon}\right) + \lambda(x^2 + y^2)^{-\alpha/2} + \gamma(x^2+y^2)^{\beta/2}.$$
The gradient of this kernel is
$$\nabla K(x,y) = 
    \begin{bmatrix}
        x\left(\dfrac{\delta\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{x^2+y^2\mathstrut}{4\varepsilon}\right) - \alpha \lambda (x^2+y^2)^{-\alpha/2-1} + \beta\gamma (x^2 + y^2)^{\beta/2-1}\right) \\
        y\left(\dfrac{\delta\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{x^2+y^2\mathstrut}{4\varepsilon}\right) - \alpha \lambda (x^2+y^2)^{-\alpha/2-1} + \beta\gamma (x^2 + y^2)^{\beta/2-1}\right)
    \end{bmatrix}.$$
If *α* = 0, we replace
*λ*(*x*<sup>2</sup>+*y*<sup>2</sup>)<sup>−*α*/2</sup> with
$$\lambda\log\left(\frac{1}{\sqrt{x^2+y^2}}\right) = -\frac{\lambda}{2}\log(x^2+y^2).$$
This gives
$$K\_{\alpha = 0}(x,y) = -\frac{\delta}{4\sqrt{\pi}\varepsilon^{3/2}}\exp\left(-\dfrac{x^2+y^2\mathstrut}{4\varepsilon}\right) - \frac{\lambda}{2}\log(x^2+y^2) + \gamma(x^2+y^2)^{\beta/2}$$
and
$$\nabla K\_{\alpha = 0}(x,y) = 
    \begin{bmatrix}
        x\left(\dfrac{\delta\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{x^2+y^2\mathstrut}{4\varepsilon}\right) - \dfrac{\lambda\mathstrut}{x^2+y^2} + \beta\gamma (x^2 + y^2)^{\beta/2-1}\right) \\
        y\left(\dfrac{\delta\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{x^2+y^2\mathstrut}{4\varepsilon}\right) - \dfrac{\lambda\mathstrut}{x^2+y^2} + \beta\gamma (x^2 + y^2)^{\beta/2-1}\right)
    \end{bmatrix}.$$
Similarly, if *β* = 0, we replace
*γ*(*x*<sup>2</sup>+*y*<sup>2</sup>)<sup>*β*/2</sup> with
$$\gamma\log\left(\sqrt{x^2+y^2}\right) = \frac{\gamma}{2}\log\left(x^2+y^2\right).$$
This gives
$$K\_{\beta = 0}(x,y) = -\frac{\delta}{4\sqrt{\pi}\varepsilon^{3/2}}\exp\left(-\dfrac{x^2+y^2\mathstrut}{4\varepsilon}\right) - \lambda\log(x^2+y^2) + \frac{\gamma}{2}\log(x^2+y^2)$$
and
$$\nabla K\_{\beta = 0}(x,y) = 
    \begin{bmatrix}
        x\left(\dfrac{\delta\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{x^2+y^2\mathstrut}{4\varepsilon}\right) - \alpha \lambda (x^2+y^2)^{-\alpha/2-1} + \dfrac{\gamma\mathstrut}{x^2+y^2}\right) \\
        y\left(\dfrac{\delta\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{x^2+y^2\mathstrut}{4\varepsilon}\right) - \alpha \lambda (x^2+y^2)^{-\alpha/2-1} + \dfrac{\gamma\mathstrut}{x^2+y^2}\right)
    \end{bmatrix}.$$

**Elliptical Kernel**

Let *a* and *b* be positive real numbers. The “elliptical kernel" is
defined as
$$K(x,y) = -\frac{\delta}{4\sqrt{\pi}\varepsilon^{3/2}}\exp\left(-\dfrac{a^2x^2+b^2y^2\mathstrut}{4\varepsilon}\right) + \lambda(x^2 + y^2)^{-\alpha/2} + \gamma(x^2+y^2)^{\beta/2}.$$
The gradient of this kernel is
$$\nabla K(x,y) = 
    \begin{bmatrix}
        x\left(\dfrac{\delta a^2\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^2x^2+b^2y^2\mathstrut}{4\varepsilon}\right) - \alpha \lambda (x^2+y^2)^{-\alpha/2-1} + \beta\gamma (x^2 + y^2)^{\beta/2-1}\right) \\
        y\left(\dfrac{\delta b^2\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^2x^2+b^2y^2\mathstrut}{4\varepsilon}\right) - \alpha \lambda (x^2+y^2)^{-\alpha/2-1} + \beta\gamma (x^2 + y^2)^{\beta/2-1}\right)
    \end{bmatrix}.$$
If *α* = 0, we replace
*λ*(*x*<sup>2</sup>+*y*<sup>2</sup>)<sup>−*α*/2</sup> with
$$\lambda\log\left(\frac{1}{\sqrt{x^2+y^2}}\right) = -\frac{\lambda}{2}\log(x^2+y^2).$$
This gives
$$K\_{\alpha = 0}(x,y) = -\frac{\displaystyle \delta}{\displaystyle 4\sqrt{\pi}\varepsilon^{3/2}}\exp\left(-\dfrac{a^2x^2+b^2y^2\mathstrut}{4\varepsilon}\right) - \frac{\lambda}{2}\log(x^2+y^2) + \gamma(x^2+y^2)^{\beta/2}$$
and
$$\nabla K\_{\alpha = 0}(x,y) = 
    \begin{bmatrix}
        x\left(\dfrac{\delta a^2\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^2x^2+b^2y^2\mathstrut}{4\varepsilon}\right) - \dfrac{\lambda\mathstrut}{x^2+y^2} + \beta\gamma (x^2 + y^2)^{\beta/2-1}\right) \\
        y\left(\dfrac{\delta b^2\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^2x^2+b^2y^2\mathstrut}{4\varepsilon}\right) - \dfrac{\lambda\mathstrut}{x^2+y^2} + \beta\gamma (x^2 + y^2)^{\beta/2-1}\right)
    \end{bmatrix}.$$
Similarly, if *β* = 0, we replace
*γ*(*x*<sup>2</sup>+*y*<sup>2</sup>)<sup>*β*/2</sup> with
$$\gamma\log\left(\sqrt{x^2+y^2}\right) = \frac{\gamma}{2}\log\left(x^2+y^2\right).$$
This gives
$$K\_{\beta = 0}(x,y) = -\frac{\delta}{4\sqrt{\pi}\varepsilon^{3/2}}\exp\left(-\dfrac{a^2x^2+b^2y^2\mathstrut}{4\varepsilon}\right) - \lambda\log(x^2+y^2) + \frac{\gamma}{2}\log(x^2+y^2)$$
and
$$\nabla K\_{\beta = 0}(x,y) = 
    \begin{bmatrix}
        x\left(\dfrac{\delta a^2\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^2x^2+b^2y^2\mathstrut}{4\varepsilon}\right) - \alpha \lambda (x^2+y^2)^{-\alpha/2-1} + \dfrac{\gamma\mathstrut}{x^2+y^2}\right) \\
        y\left(\dfrac{\delta b^2\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^2x^2+b^2y^2\mathstrut}{4\varepsilon}\right) - \alpha \lambda (x^2+y^2)^{-\alpha/2-1} + \dfrac{\gamma\mathstrut}{x^2+y^2}\right)
    \end{bmatrix}.$$

**Elliptical *c*-Norm Kernel**

Let *c* be an even positive integer. The “elliptical *c*-norm kernel" is
defined as
$$K(x,y) = -\frac{\delta}{4\sqrt{\pi}\varepsilon^{3/2}}\exp\left(-\dfrac{a^cx^c+b^cy^c\mathstrut}{4\varepsilon}\right) + \lambda(x^c + y^c)^{-\alpha/c} + \gamma(x^c+y^c)^{\beta/c}.$$
The gradient of this kernel is
$$\nabla K(x,y) = \begin{bmatrix}
        x^{c-1}\left(\dfrac{\delta c a^c\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^cx^c+b^cy^c\mathstrut}{4\varepsilon}\right)-\alpha\lambda(x^c+y^c)^{-\alpha/c-1} + \gamma\beta(x^c+y^c)^{\beta/c-1}\right) \\
        y^{c-1}\left(\dfrac{\delta c b^c\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^cx^c+b^cy^c\mathstrut}{4\varepsilon}\right)-\alpha\lambda(x^c+y^c)^{-\alpha/c-1} + \gamma\beta(x^c+y^c)^{\beta/c-1}\right)
    \end{bmatrix}.$$
If *α* = 0, we replace
*λ*(*x*<sup>*c*</sup>+*y*<sup>*c*</sup>)<sup>−*α*/*c*</sup> with
$$\lambda\log\left((x^c+y^c)^{-1/c}\right) = -\frac{\lambda}{c}\log(x^c+y^c).$$
This gives
$$K\_{\alpha = 0}(x,y) = -\frac{\delta}{4\sqrt{\pi}\varepsilon^{3/2}}\exp\left(-\dfrac{a^cx^c+b^cy^c\mathstrut}{4\varepsilon}\right)  -\frac{\lambda}{c}\log(x^c+y^c) + \gamma(x^c+y^c)^{\beta/c}$$
and
$$\nabla K\_{\alpha = 0}(x,y) = \begin{bmatrix}
        x^{c-1}\left(\dfrac{\delta c a^c\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^cx^c+b^cy^c\mathstrut}{4\varepsilon}\right)-\dfrac{\lambda\mathstrut}{x^c+y^c} + \gamma\beta(x^c+y^c)^{\beta/c-1}\right) \\
        y^{c-1}\left(\dfrac{\delta c b^c\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^cx^c+b^cy^c\mathstrut}{4\varepsilon}\right)-\dfrac{\lambda\mathstrut}{x^c+y^c} + \gamma\beta(x^c+y^c)^{\beta/c-1}\right)
    \end{bmatrix}.$$
If *β* = 0, we replace
*γ*(*x*<sup>*c*</sup>+*y*<sup>*c*</sup>)<sup>*β*/*c*</sup> with
$$\gamma\log\left((x^c+y^c)^{1/c}\right) = \frac{\gamma}{c}\log(x^c+y^c).$$
This gives
$$K\_{\beta = 0}(x,y) = -\frac{\delta}{4\sqrt{\pi}\varepsilon^{3/2}}\exp\left(-\dfrac{a^cx^c+b^cy^c\mathstrut}{4\varepsilon}\right) + \lambda(x^c + y^c)^{-\alpha/c} + \frac{\gamma}{c}\log(x^c+y^c)$$
and
$$\nabla K\_{\beta = 0}(x,y) = \begin{bmatrix}
        x^{c-1}\left(\dfrac{\delta c a^c\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^cx^c+b^cy^c\mathstrut}{4\varepsilon}\right)-\alpha\lambda(x^c+y^c)^{-\alpha/c-1} + \dfrac{\gamma\mathstrut}{x^c+y^c}\right) \\
        y^{c-1}\left(\dfrac{\delta c b^c\mathstrut}{8\sqrt{\pi}\varepsilon^{5/2}}\exp\left(-\dfrac{a^cx^c+b^cy^c\mathstrut}{4\varepsilon}\right)-\alpha\lambda(x^c+y^c)^{-\alpha/c-1} + \dfrac{\gamma\mathstrut}{x^c+y^c}\right)
    \end{bmatrix}.$$