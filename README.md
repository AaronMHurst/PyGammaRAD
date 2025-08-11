# PyGammaRAD

This project is a Python implementation of a library that can be used in the calculation of &gamma;-ray angular distribution coefficients in addition to a general purpose angular momentum calculator for the evaluation of quantities that underly the determination of said coefficients.  Stritcly speaking, it is only Clebsch-Gordan and Racah coefficients that are needed in the theoretical descrition of the angular distribution functions considered here.  However, given to the close relationship these coefficients have with other angular mementum symbols typically used to describe coupling and recoupling schemes in quantum mechanical applications involving angular momenta, to make the package more complete and broaden its utility we also provide methods to readily evaluate the Wigner 3-j, 6-j, and 9-j symbols.  Finally, it is also intended that this software package can serve as an API to the methods and data tables descibed in the *"Tables of Coefficients for Angular Distribution of Gamma Rays from Aligned Nuclei"* [[1]](#1).

![W 56Fe](W_functions_56Fe.png?raw=True "Angular distribution functions for three different transitions in <sup>56</sup>Fe")

## References

<a id="1">[1]</a>
T. Yamazaki,
*"Tables of Coefficients for Angular Distribution of Gamma Rays from Aligned Nuclei"*,
Nucl. Data, Sect. A, Vol. **3**, Num. 1 (1967).
