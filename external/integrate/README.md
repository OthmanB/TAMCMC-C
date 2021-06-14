# integrate Glm~F(theta, phi).|Ylm|^2
Small program for integrating Glm(theta, phi), term that can be used to define activity effect on mode splittings. It uses modified functions taken from https://github.com/CD3/libIntegrate to perform integration in C++ using Gaussian Legendre Quadrature

This code can calculate Glm as defined in Gizon 2002, AN, 323, 251–253
But can be used to have more flexible activity zone hypothesis than that paper.

# How to use
Use cmake to compile into a build directory:
       'mkdir build'
       'cd build'
       'cmake ..'
       'make'

This will generate an executable with parameters as per defined into the main() of ylm.cpp
Note that this main is just for tests. A proper usage would require to include ylm.cpp into a larger projet,
after commenting the main().


PS: This version has the main() section of activtiy.cpp commented because it is intended to be used as a plug-in for the TAMCMC code.
If you want to run it as a standalone, uncomment main() section of your choice in activity.cpp.

Othman Benomar
