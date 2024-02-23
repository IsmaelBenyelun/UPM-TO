# UPM-TO

Julia package to run the "Updated Properties Model" (Saucedo et al., 2023; Ben-Yelun et al. 2023) [1, 2] topology optimization algorithm.

To load the corresponding (local) packages `tfgfem` and `optim`, open `julia` in this directory and execute `dev tfgfem` and `dev optim` in the Julia package manager (pressing `]` in the REPL)

```
julia> ]
(@vX.X) pkg> dev tfgfem

(@vX.X) pkg> dev optim
```

[1] Saucedo-Mora, L., Ben-Yelun, I., García-Modet, H., Sanz-Gómez, M. Á., & Montans, F. J. (2023). The Updated Properties Model (UPM): A topology optimization algorithm for the creation of macro–micro optimized structures with variable stiffness. Finite Elements in Analysis and Design, 223, 103970.

[2] Ben-Yelun, I., Saucedo-Mora, L., Sanz, M. Á., Benítez, J. M., & Montans, F. J. (2023). Topology optimization approach for functionally graded metamaterial components based on homogenization of mechanical variables. Computers & Structures, 289, 107151.
