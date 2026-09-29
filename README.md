# Numerical Analysis Coursework

Jupyter notebooks for numerical analysis coursework, covering interpolation, numerical differentiation and integration, and the numerical solution of ordinary differential equations. Notebooks are organized by course week and include CSV data files where needed.

## Contents

### Week 6: Interpolation

- Equipment: [Chebyshev Nodes](Week6-Interpolation/pp0_Equipment/pa1_ChebyshevNode.ipynb), [Horner Table](Week6-Interpolation/pp0_Equipment/pa2_HornerTable.ipynb)
- [Lagrange Interpolation](Week6-Interpolation/pp1_Lagrange/pp1_LagrangeInterpolation.ipynb)
- Newton: [Interpolation](Week6-Interpolation/pp2_Newton/pp2a_NewtonInterpolation.ipynb), [Fixed Gap](Week6-Interpolation/pp2_Newton/pp2b_NewtonFixedGap.ipynb)
- Centered differences: [Gauss I](Week6-Interpolation/pp3_CenteredDifference/pp3a_GaussI.ipynb), [Gauss II](Week6-Interpolation/pp3_CenteredDifference/pp3b_GaussII.ipynb), [Stirling](Week6-Interpolation/pp3_CenteredDifference/pp3c_Stirling.ipynb), [Bessel](Week6-Interpolation/pp3_CenteredDifference/pp3d_Bessel.ipynb)
- Inverse interpolation: [Inverse Function](Week6-Interpolation/pp4_InverseInterpolation/pp4a_InverseFunction.ipynb), [Loop Inverse](Week6-Interpolation/pp4_InverseInterpolation/pp4b_LoopInverse.ipynb)
- Splines: [Linear](Week6-Interpolation/pp5_SplineInterpolation/pp5a_LinearSpline.ipynb), [Quadratic](Week6-Interpolation/pp5_SplineInterpolation/pp5b_QuadraticSpline.ipynb), [Cubic](Week6-Interpolation/pp5_SplineInterpolation/pp5c_CubicSpline.ipynb)
- [Least Sum of Squares](Week6-Interpolation/pp6-LeastSumSquare/pp6_LSQ.ipynb)

### Week 7: Differentiation and Integration

- [Lagrange Differentiation](Week7-Integral/pp1_Derivative/pp1_LagrangeDerivative.ipynb)
- Definite integrals: [Lagrange](Week7-Integral/pp2_DefiniteIntegral/pp2a_LagrangeIntegral.ipynb), [Trapezoidal Rule](Week7-Integral/pp2_DefiniteIntegral/pp2b_Trapezoid.ipynb), [Simpson's Rule](Week7-Integral/pp2_DefiniteIntegral/pp2c_Simpson.ipynb), [Newton-Cotes](Week7-Integral/pp2_DefiniteIntegral/pp2d_NewtonCotez.ipynb)

### Week 8: Ordinary Differential Equations

- Euler methods: [Forward](Week8-Differential/pp1a_EulerForward.ipynb), [Backward](Week8-Differential/pp1b_EulerBackward.ipynb)
- [Modified Heun](Week8-Differential/pp1c_HeunModified.ipynb)
- [Runge-Kutta](Week8-Differential/pp2_RungeKetta.ipynb)
- Multistep methods: [Adams-Bashforth Explicit](Week8-Differential/pp3a_ABExplicit.ipynb), [Adams-Moulton Implicit](Week8-Differential/pp3b_AMImplicit.ipynb)

## Requirements

- Python 3
- JupyterLab or Jupyter Notebook
- NumPy, pandas, Matplotlib, and SymPy

Install the required packages with:

```bash
python -m pip install jupyter numpy pandas matplotlib sympy
```

From the repository directory, start JupyterLab with:

```bash
jupyter lab
```

Then open a notebook from the relevant week folder. Keep each notebook alongside its referenced CSV data files when running it locally.

## Google Colab

[Open the project notebook in Google Colab](https://colab.research.google.com/drive/1CHWoonrCRVfp68gAgrNehbeefCuciOEg?usp=sharing).

