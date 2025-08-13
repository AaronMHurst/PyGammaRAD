# PyGammaRAD

The `PyGammaRAD` project (*Python project for Gamma-Ray Angular Distributions*) is a Python implementation of a library that can be used in the calculation of &gamma;-ray angular distribution coefficients in addition to a general purpose angular momentum calculator for the evaluation of quantities that underly the determination of said coefficients.  Stritcly speaking, it is only Clebsch-Gordan and Racah coefficients that are needed in the theoretical description of the angular distribution functions considered here.  However, given the close relationship these coefficients share with other angular momentum symbols typically used to describe coupling and recoupling schemes in quantum mechanical applications involving angular momenta, to make the package more complete and broaden its utility we also provide methods to readily evaluate the Wigner *3-j*, *6-j*, and *9-j* symbols.  Finally, it is also intended that this software package can serve as an API to the methods and nuclear data tables descibed in the *"Tables of Coefficients for Angular Distribution of Gamma Rays from Aligned Nuclei"* [[1]](#1).  The figure below shows an example of the overall angular distribution functions used to describe three different &gamma;-ray transitions in <sup>56</sup>Fe assuming complete nuclear alignment.

![W 56Fe](W_functions_56Fe.png?raw=True "Angular distribution functions for three different transitions in <sup>56</sup>Fe")

Here, the angular distribution function is written as

$$ W(\theta) = 1 + A_{2}P_{2}(\cos\theta) +A_{4}P_{4}(\cos\theta), $$

where *A<sub>k</sub>* are the theoretical anisotropy coefficients for transitions in aligned nuclei and *P<sub>k</sub>* are the corresponding Legendre polynomials of a given order *k*.

## Building and installation

After cloning the repository, this project can then be built and installed by running the `installation.sh` script at the terminal command line of the project directory:

```Bash
$ git clone https://github.com/AaronMHurst/PyGammaRAD.git
$ cd pace_ensdf
$ sh installation.sh
```

In the future, the project will also be deployed to the PyPI repository.


## Testing

A suite of Python modules containing 123 unit tests have been written for this project.  These unit-test scripts are located in the `tests` folder.  To run the test suite and ensure they work with the local Python environment, run `tox` at the command line of the project directory containing the `tox.ini` file:

```Bash
$ tox -r
```

This project has the following Python-package dependencies: `numpy`, `pandas`, and `pytest`.  The test session is automatically started after building against the required virtual Python environment.


## Running the software

Subsequent to the installation process, the `PyGammaRAD` libray is accessibles from any location by importing the package and creating an instance of the `AngularMomentum` class:

```Bash
$ python
```
```python
>>> import PyGammaRAD as pg
>>> am = pg.AngularMomentum()
```

To help illustrate the workflow and utility of the software, inclduing verification of methods againts published results and tabulated data, the project comes with three different `Jupyter Notebooks` for the user to run through:

* `am_coupling_check.ipynb`: Provides a useful guide for executing angular momentum methods involving Clebsch-Gordan and Racah coefficients in addition to the Wigner *3-j*, *6-j*, and *9-j* symbols.  Well-known published results are compared and verified in this `Notebook`.

* `angular_distributions`: This `Notebook` illustrates the use of the methods involved in the calculation of the overall &gamma;-ray angular distribution function including Legendre polynomials available to the package.  For comparison, experimental anisotropy-attenuation coefficients for transitions in <sup>56</sup>Fe are compared to the theoretically-deduced results.

* `yamazaki_tables`: This `Notebook` serves as a check of the original angular distribution tensors and coefficients published and tabulated by Yamazaki [[1]](#1).


## Docstrings

All `PyGammaRAD` classes and functions have supporting docstrings.  Please refer to the individual dosctrings for more information on any particular function including how to use it.  The dosctrings for each method generally have the following structure:

* A short explanation of the function.
* A list and description of arguments that need to be passed to the function.
* The return value of the function.
* An example(s) invoking use of the function.

To retrieve a list of all available methods simply execute the following command in a Python interpreter:

```python
>>> help(am)
```

Or, to retrieve the docstring for a particular method, e.g., the callable `symb6j` to evaluate the corresponding Wigner *6-j* symbol:

```python
>>> help(am.symb6j)
```

The `Jupyter Notebooks` also illustrate docstring retrieval for certain methods.

## Summary of angular distribution functions and methods

The table below summarizes the angular distribution functions given in the reference article by Yamazaki [[1]](#1) and their corresponding callable methods available to the `PyGammaRAD` software package.  The relevant arguments, listed in order where needed, are defined as:

* *k* : Order of the coefficient or polynomial degree.
* *J<sub>i</sub>* : Initial nuclear level of the associated &gamma;-ray transition.
* *J<sub>f</sub>* : Final nuclear level of the associated &gamma;-ray transition.
* *L<sub>1</sub>* : First multipole order.
* *L<sub>2</sub>* : Second multipole order.
* &delta;<sub>&gamma;</sub> : &gamma;-ray multipole mixing ratio.
* *A<sub>k</sub>* : Anisotropy coefficient of order *k*.

| Function [[1]](#1) | Method | Arguments |
| --- | --- | --- |
| Equation (2) | `dist_W` | *A<sub>k</sub>* |
| Equation (4) | `calc_F` | *k*, *J<sub>f</sub>*, *L<sub>1</sub>*, *L<sub>2</sub>*, *J<sub>i</sub>* |
| Equation (6) | `calc_B` | *k*, *J* |
| Equation (7) | `A_max` | *k*, *J<sub>i</sub>*, *L<sub>1</sub>*, *L<sub>2</sub>*, *J<sub>f</sub>*, &delta;<sub>&gamma;</sub> |
| Equation (8) | `calc_BF` | *k*, *J<sub>f</sub>*, *L<sub>1</sub>*, *L<sub>2</sub>*, *J<sub>i</sub>* |
| Equation (13) | `U_coeff` | *k*, *J<sub>i</sub>*, *L<sub>1</sub>*, *L<sub>2</sub>*, *J<sub>f</sub>*, &delta;<sub>&gamma;</sub> |
| Equation (14) | `calc_u` | *k*, *J<sub>i</sub>*, *L<sub>1</sub>*, *J<sub>f</sub>* |

## Summary angular momentum methods

The set of angular momentum functions and callable methods available to `PyGammaRAD` is tabulated below.  The required arguments are listed in the order in which they should be passed to their corresponding method.

| Name | Coefficient/Symbol | Method | Arguments |
| --- | --- | --- | --- |
|Clebsch-Gordan | $<j_{1} m_{1} j_{2} m_{2} \|j m>$ | `cg` | $j_{1}$, $m_{1}$, $j_{2}$, $m_{2}$, $j$, $m$ |
| Wigner 3-*j* | $( j_{1} \quad j_{2} \quad j$ <br> $m_{1} \quad m_{2} \quad m )$  | `symb3j` | $j_{1}$, $j_{2}$, $j$, $m_{1}$, $m_{2}$, $m$ |
| Racah | $W(j_{1} j_{2} j_{3} j_{4}; j_{5} j_{6})$ | `racah` | $j_{1}$, $j_{2}$, $j_{3}$, $j_{4}$, $j_{5}$, $j_{6}$ |
| Wigner 6-*j* | ${ j_{1} \quad j_{2} \quad j_{3}$ <br> $j_{4} \quad j_{5} \quad j_{6} }$  | `symb3j` | $j_{1}$, $j_{2}$, $j_{3}$, $j_{4}$, $j_{5}$, $j_{6}$ |


## References

<a id="1">[1]</a>
T. Yamazaki,
*"Tables of Coefficients for Angular Distribution of Gamma Rays from Aligned Nuclei"*,
Nucl. Data, Sect. A, Vol. **3**, Num. 1 (1967).
