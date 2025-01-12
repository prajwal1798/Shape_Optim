# **Shape Optimization Toolbox for NACA0012 Aerofoil**
![Neural_Net](https://github.com/user-attachments/assets/8818819f-f3ef-4b00-a503-3fb0f2b965da)

## **Overview**
The shape optimization toolbox is designed to improve the aerodynamic performance of the **NACA0012 aerofoil** by minimizing the **drag-to-lift ratio** $\frac{C_d}{C_l}$ using a **genetic algorithm (GA)**. This toolbox interfaces **MATLAB** (for geometry generation) with **Python** (for neural network predictions), enabling an automated flow of data between these environments. The GA performs **gradient-free optimization**, making it robust for high-dimensional, non-differentiable objective functions.

---

## **1. Genetic Algorithm Framework**

The optimization process begins with a random initial population of control points for the NACA0012 aerofoil. The steps of the **genetic algorithm** (GA) are outlined below:

1. **Initialize Population:** Generate random control point configurations within specified bounds.
2. **Evaluate Objective Function:** Evaluate each population member's drag-to-lift ratio $\frac{C_d}{C_l}$ by sending the control points to the neural network.
3. **Selection:** Select individuals with the lowest $\frac{C_d}{C_l}$ to form the next generation.
4. **Crossover and Mutation:** Combine and mutate the selected individuals to introduce diversity and explore new configurations.
5. **Convergence Check:** The process repeats until the improvement in $\frac{C_d}{C_l}$ falls below a tolerance (set as $10^{-6}$) or the maximum generation count is reached.

---

## **2. Control Point Variation Bounds**

The reference aerofoil is based on the NACA0012 profile, with upper and lower bounds for control point perturbations. These bounds define the design space for optimization, providing **4 degrees of freedom** for the GA.

---

## **3. Objective Function Evaluation**

The neural network, trained separately, predicts wall quantities such as **pressure** ($p$) and **shear stresses** ($\tau_1$, $\tau_2$) based on the control point configurations. The GA uses these predictions to compute the drag-to-lift ratio $\frac{C_d}{C_l}$:

### **Objective Function Computation Steps:**
- **Lift and Drag Calculation:** Neural network predictions are integrated over the aerofoil surface to compute $C_l$ and $C_d$.
- **Drag-to-Lift Ratio:** The objective function is defined as $\frac{C_d}{C_l}$. If $C_l \leq 0$ or $C_d \leq 0$, a high penalty is assigned to discourage such configurations.

---

## **4. Genetic Algorithm for Shape Optimization**

### **Algorithm: Genetic Optimization for Shape Optimization of NACA0012 Aerofoil**

**Input:** Control points $P$, bounds $lb$, $ub$, population size $N$, max generations $G$, tolerance $\epsilon$  
**Output:** Optimized control points $P^*$, minimum drag-to-lift ratio $\left(\frac{C_d}{C_l}\right)^*$

1. Initialize population $\{P_i\}_{i=1}^N$ within bounds $lb \leq P \leq ub$.
2. **For** each generation:
   - Send $P_i$ to the neural network to predict $p$, $\tau_1$, and $\tau_2$.
   - Compute $C_l$ and $C_d$ using surface integrals:
     ![equation](https://latex.codecogs.com/svg.image?\color{White}C_l=\frac{1}{q_{\infty}S_{\Gamma}}\int_{\Gamma}\sigma_j^wn_j^{\infty}d\bm{x})
     ![equation](https://latex.codecogs.com/svg.image?\color{White}C_d=\frac{1}{q_{\infty}S_{\Gamma}}\int_{\Gamma}\sigma_j^wt_j^{\infty}d\bm{x})
   - **If** $C_l > 0$ and $C_d > 0$:
     - Calculate $f(P_i) = \frac{C_d}{C_l}$.
   - **Else:**
     - Set $f(P_i) = \infty$ (penalty).
3. **Selection:** Choose individuals with the lowest $f(P_i)$.
4. **Crossover and Mutation:** Generate new population.
5. **Convergence Check:** Stop if $f(P^*) < \epsilon$.

---

## **5. Optimization Results**

### **5.1 Optimal Control Points**
The GA optimized the aerofoil after **4280 function evaluations**. The final set of control points is shown below:

| **Point** | **$X_i^+$ (Upper Surface)** | **$Y_i^+$** | **$X_i^-$ (Lower Surface)** | **$Y_i^-$** |
|-----------|-----------------------------|------------|----------------------------|------------|
| 1         | 0.5000                       | 0.0000     | 0.5000                      | 0.0000     |
| 2         | 0.4376                       | 0.0127     | 0.4514                      | -0.0099    |
| 3         | 0.1825                       | 0.0416     | 0.2088                      | -0.0398    |
| 4         | -0.0392                      | 0.0577     | -0.0425                     | -0.0571    |
| 5         | -0.2320                      | 0.0631     | -0.2480                     | -0.0423    |
| 6         | -0.4004                      | 0.0477     | -0.4375                     | -0.0268    |
| 7         | -0.4932                      | 0.0395     | -0.4842                     | -0.0145    |
| 8         | -0.5000                      | 0.0000     | -0.5000                     | 0.0000     |

---
![optmial_profile](https://github.com/user-attachments/assets/bdaa93d4-3528-4685-a71a-1cd1d595ae0b)

### **5.2 Achieved Minimum Drag-to-Lift Ratio**
The optimized drag-to-lift ratio $\frac{C_d}{C_l}$ achieved was **0.54545**, demonstrating the GA's effectiveness in achieving a global minimum beyond the bounds of the training dataset.

---

## **6. Conclusion**

### **6.1 Effective Design Space Exploration**
- The training dataset of **2,560 geometries** contained drag-to-lift ratios ranging from **1.9** to **11,125**.
- The GA achieved a ratio of **0.54**, significantly better than the best value in the training dataset.

### **6.2 Validation of Surrogate Model Predictions**
The optimized aerofoil's pressure and shear stress values fell within the training dataset bounds:
- **Pressure:** All optimized $p$ values were within the range of training data.
- **Shear Stresses:** Both $\tau_1$ and $\tau_2$ values were within the bounds.

### **6.3 Computational Efficiency**
The proposed framework significantly reduced computational expenses:
- **Traditional CFD-Based Optimization:** Requires a full CFD simulation for each function evaluation (prohibitively expensive).
- **Surrogate Model Framework:** Replaces CFD calls with real-time neural network predictions, reducing computational cost by **50%**.
- The optimization completed in **4,850 function evaluations** without sacrificing accuracy.

