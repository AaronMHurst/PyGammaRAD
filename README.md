# PyGammaRAD

The `PyGammaRAD` project (*Python project for Gamma-Ray Angular Distributions*) is a Python implementation of a library that can be used in the calculation of &gamma;-ray angular distribution coefficients in addition to a general purpose angular momentum calculator for the evaluation of quantities that underly the determination of said coefficients.  Stritcly speaking, it is only Clebsch-Gordan and Racah coefficients that are needed in the theoretical description of the angular distribution functions considered here.  However, given the close relationship these coefficients share with other angular momentum symbols typically used to describe coupling and recoupling schemes in quantum mechanical applications involving angular momenta, to make the package more complete and broaden its utility we also provide methods to readily evaluate the Wigner *3-j*, *6-j*, and *9-j* symbols.  Finally, it is also intended that this software package can serve as an API to the methods and nuclear data tables descibed in the *"Tables of Coefficients for Angular Distribution of Gamma Rays from Aligned Nuclei"* [[1]](#1) and those in the appendix of *"Angular Distributions of Gamma Rays in Terms of Phase-Defined Reduced Matrix Elements"* [[2]](#2).  The figure below shows an example of the overall angular distribution functions used to describe three different &gamma;-ray transitions in <sup>56</sup>Fe assuming complete nuclear alignment.

![W 56Fe](W_functions_56Fe.png?raw=True "Angular distribution functions for three different transitions in <sup>56</sup>Fe")

Here, the angular distribution function is written as

$$ W(\theta) = 1 + A_{2}P_{2}(\cos\theta) +A_{4}P_{4}(\cos\theta), $$

where *A<sub>k</sub>* are the theoretical anisotropy coefficients for transitions in aligned nuclei and *P<sub>k</sub>* are the corresponding Legendre polynomials of a given order *k*.

## Building and installation

After cloning the repository, this project can then be built and installed by running the `installation.sh` script at the terminal command line of the project directory:

```Bash
$ git clone https://github.com/AaronMHurst/PyGammaRAD.git
$ cd PyGammaRAD
$ sh installation.sh
```

In the future, the project will also be deployed to the PyPI repository.


## Testing

A suite of Python modules containing 140 unit tests have been written for this project.  These unit-test scripts are located in the `tests` folder.  To run the test suite and ensure they work with the local Python environment, run `tox` at the command line of the project directory containing the `tox.ini` file:

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

To help illustrate the workflow and utility of the software, inclduing verification of methods against published results and tabulated data, the project comes with four different `Jupyter Notebooks` for the user to run through:

* `am_coupling_check.ipynb`: Provides a useful guide for executing angular momentum methods involving Clebsch-Gordan and Racah coefficients in addition to the Wigner *3-j*, *6-j*, and *9-j* symbols.  Well-known published results are compared and verified in this `Notebook`.

* `angular_distributions.ipynb`: This `Notebook` illustrates the use of the methods involved in the calculation of the overall &gamma;-ray angular distribution function including Legendre polynomials available to the package.  For comparison, experimental anisotropy-attenuation coefficients for transitions in <sup>56</sup>Fe are also compared to the theoretically-deduced results.

* `yamazaki_tables.ipynb`: This `Notebook` serves as a check of the original angular distribution tensors and coefficients published and tabulated by Yamazaki [[1]](#1).

* `rose_brink_tables.ipynb`: This `Notebook` serves as a check of the original angular distribution tensors and coefficients published and tabulated by Rose and Brink [[2]](#2).


## Docstrings

All `PyGammaRAD` classes and functions have supporting docstrings.  Please refer to the individual dosctrings for more information on any particular function including how to use it.  The dosctrings for each method generally have the following structure:

* A short explanation of the function.
* A list and description of arguments that need to be passed to the function.
* The return value of the function.
* Exceptions that may be raised.
* An example(s) invoking use of the function.

To retrieve a list of all available methods simply execute the following command in a Python interpreter:

```python
>>> help(am)
```

Or, to retrieve the docstring for a particular method, e.g., the callable `symb6j` to evaluate the corresponding Wigner *6-j* symbol:

```python
>>> help(am.symb6j)
```

The `Jupyter Notebooks` provided also illustrate docstring retrieval for certain methods.

## Summary of angular distribution functions and methods

The table below summarizes the angular distribution functions given in the reference articles by Yamazaki [[1]](#1) and Rose and Brink [[2]](#2) and their corresponding callable methods available within the `PyGammaRAD` software package.  The relevant arguments, listed in order where needed, are defined as:

* *k* : Order of the coefficient or polynomial degree.
* *J<sub>i</sub>* : Initial nuclear level of the associated &gamma;-ray transition.
* *J<sub>f</sub>* : Final nuclear level of the associated &gamma;-ray transition.
* *L<sub>1</sub>* : First multipole order.
* *L<sub>2</sub>* : Second (interfering) multipole order.
* &delta;<sub>&gamma;</sub> : &gamma;-ray multipole mixing ratio.
* *A<sub>k</sub>* : Anisotropy coefficient of order *k*.
* *l1*: Orbital angular momentum of first partial wave.
* *l2*: Orbital angular momentum of second (interfering) partial wave.
* *s*: Reaction channel spin.
* *m*: Magnetic substate projection quantum number.

Optional arguments are preceded by an asterisk (<sup>*</sup>).  Note that in the case of a mixed transition *L<sub>1</sub> &ne; L<sub>2</sub>*, while in the case of a pure stretched transition *L<sub>1</sub> = L<sub>2</sub>*.  For a given transition the interfering multipole order is defined such that *L<sub>2</sub> = L<sub>1</sub> + 1*.

| Returned quantity | Function [[1]](#1),[[2]](#2) | Method | Arguments |
| --- | --- | --- | --- |
| $W(\theta)$ | Equation (2) [[1]](#1) | `dist_W` | *A<sub>k</sub>*, <sup>*</sup> *[&theta;]* |
| $F_{k}(J_{f}L_{1}L_{2}J_{i})$ | Equation (4) [[1]](#1) | `calc_F` | *k*, *J<sub>f</sub>*, *L<sub>1</sub>*, *L<sub>2</sub>*, *J<sub>i</sub>* |
| $B_{k}(J)$  | Equation (6) [[1]](#1) | `calc_B` | *k*, *J* |
| $A_{k}^{max}(J_{i}L_{1}L_{2}J_{f})$  | Equation (7) [[1]](#1) | `A_max` | *k*, *J<sub>i</sub>*, *L<sub>1</sub>*, *L<sub>2</sub>*, *J<sub>f</sub>*, &delta;<sub>&gamma;</sub> |
| $B_{k}(J_{i})F_{k}(J_{f}L_{1}L_{2}J_{i})$ | Equation (8) [[1]](#1) | `calc_BF` | *k*, *J<sub>f</sub>*, *L<sub>1</sub>*, *L<sub>2</sub>*, *J<sub>i</sub>* |
| $U_{k}(J_{i}L_{1}L_{2}J_{f})$ | Equation (13) [[1]](#1) | `U_coeff` | *k*, *J<sub>i</sub>*, *L<sub>1</sub>*, *L<sub>2</sub>*, *J<sub>f</sub>*, &delta;<sub>&gamma;</sub> |
| $u_{k}(J_{i}L_{1}J_{f})$  | Equation (14) [[1]](#1); Equation (3.45) [[2]](#2) | `calc_u` | *k*, *J<sub>i</sub>*, *L<sub>1</sub>*, *J<sub>f</sub>* |
| $R_{k}(L_{1}L_{2}J_{i}J_{f})$ | Equations (3.36) & (3.37) [[2]](#2) | `calc_R` | *k*, *L<sub>1</sub>*, *L<sub>2</sub>*, *J<sub>i</sub>*, *J<sub>f</sub>* |
| $S_{k}(l_{1} l_{2} J s)$ | Equation (3.59) [[2]](#2) | `calc_S` | *l<sub>1</sub>*, *l<sub>2</sub>*, *J*, *s*, <sup>*</sup> *k* |
| $\rho_{k}(J m)$ | Equation (3.63) [[2]](#2) | `calc_p` | *k*, *J*, <sup>*</sup> *m* |

*[&theta;]*: A user defined range may be passed to the $W(\theta)$ function as an optional argument; by default a range of (0<sup>o</sup>,180<sup>o</sup>) is assumed.  Although it is only the non-even $P_{0}(\cos\theta)$, $P_{2}(\cos\theta)$, and $P_{4}(\cos\theta)$ Legendre polynomials that are needed to describe the angular distribution function, to enhance the utility of the `PyGammaRAD` methods the first 11 polynomials of the Legendre series are available in the library.  These methods may be called as `lpoly<k>` for the corresponding Legendre polynomial $P_{k}$, where $k$ ranges in value from $k=0$ to $k=10$.  Again, by default an angular range of (0<sup>o</sup>,180<sup>o</sup>) is adopted for the Legendre polynomial methods although a user-defined range can also be provided as an optional argument to override the default range or its granularity.

## Summary of angular momentum methods

The set of angular momentum functions and callable methods available to the `PyGammaRAD` library is tabulated below.  The required arguments are listed in the order in which they should be passed to their corresponding method.  The arguments are defined as:

* *j* : Angular momentum vector.
* *m* : magnetic substate quantum number (i.e., *z*-axis projection).

| Returned quantity | Coefficient/Symbol | Method | Arguments |
| --- | --- | --- | --- |
|Clebsch-Gordan | $<j_{1} m_{1} j_{2} m_{2} \|j m>$ | `cg` | *j<sub>1</sub>*, *m<sub>1</sub>*, *j<sub>2</sub>*, *m<sub>2</sub>*, *j*, *m* |
| Wigner 3-*j* | $j_{1}$ $j_{2}$ $j$ <br> $m_{1}$ $m_{2}$ $m$  | `symb3j` | *j<sub>1</sub>*, *j<sub>2</sub>*, *j*, *m<sub>1</sub>*, *m<sub>2</sub>*, *m* |
| Racah | $W(j_{1} j_{2} j_{3} j_{4}; j_{5} j_{6})$ | `racah` | *j<sub>1</sub>*, *j<sub>2</sub>*, *j<sub>3</sub>*, *j<sub>4</sub>*, *j<sub>5</sub>*, *j<sub>6</sub>* |
| Wigner 6-*j* | $j_{1}$ $j_{2}$ $j_{3}$ <br> $j_{4}$ $j_{5}$ $j_{6}$ | `symb6j` | *j<sub>1</sub>*, *j<sub>2</sub>*, *j<sub>3</sub>*, *j<sub>4</sub>*, *j<sub>5</sub>*, *j<sub>6</sub>* |
| Wigner 9-*j* | $j_{1}$ $j_{2}$ $j_{3}$ <br> $j_{4}$ $j_{5}$ $j_{6}$ <br> $j_{7}$ $j_{8}$ $j_{9}$ | `symb9j` | *j<sub>1</sub>*, *j<sub>2</sub>*, *j<sub>3</sub>*, *j<sub>4</sub>*, *j<sub>5</sub>*, *j<sub>6</sub>*, *j<sub>7</sub>*, *j<sub>8</sub>*, *j<sub>9</sub>* |

Note that the ordering of *j<sub>i</sub>* in the second column **does not** reflect the relationship between the Racah coefficient and its corresponding Wigner 6-*j* symbol; the labeling of the *j<sub>i</sub>* vectors is intended to facilitate ease of user input with the passing of arguments from *left-to-right*, *top-to-bottom*, as they appear in the written form of the coefficient/symbol.

The following conditions apply when handling Clebsch-Gordan coefficients and Wigner 3-*j* symbols:

* All *m* quantum numbers must be projections of their respective *j* values, i.e., for integral *j*, all *m* projections must also be integral, and likewise, half-integral *j* must also have corresponding half-integral *m* projections.  Mathematically, this rule can be expressed as
  * If 2*j<sub>i</sub>* $\bmod$ 2 = 0, then 2*m<sub>i</sub>* $\bmod$ 2 = 0;
  * If 2*j<sub>i</sub>* $\bmod$ 2 = 1, then 2*m<sub>i</sub>* $\bmod$ 2 = 1.
* Each *m* projection must satisfy the relation |*m<sub>i</sub>*| $\leq$ *j<sub>i</sub>*.
* *m<sub>1</sub>* + *m<sub>2</sub>* = m.


Additionally, for all coefficents and symbols listed above, the angular momentum vectors must satisfy the triangle inequality theorem in order to form a (*j<sub>1</sub>*, *j<sub>2</sub>*, *j<sub>3</sub>*) triad arising from anti-parallel and parallel coupling mechanisms:

* |*j<sub>1</sub>* - *j<sub>2</sub>*| $\leq$ *j<sub>3</sub>* $\leq$ *j<sub>1</sub>* + *j<sub>2</sub>*.

## Summary of Table API methods

The following set of methods enable user retrieval and manipulation of the data presented in Table 1, Table 2(a), and Table 2(b) of the original work by Yamazaki [[1]](#1).  The arguments, where required, are again listed in the order in which they should passed to their respective methods.  All physical quantities have their usual meanings defined earlier.  A few notes regarding the optional arguments and limitations on other certain arguments:

* `get_table1` : Method may be called (i) without any arguments to return both integral-*J* and half-integral *J* results, (ii) by passing `0` to return integral-*J* results only, or (iii) by passing `0.5` to return half-integral *J* results only.
* `get_row_table2` : Method only takes values of `2` or `4` as integer arguments for *k*; acceptable key-word arguments are `coeff='F'` cf. Equation (4) [[1]](#1), `coeff='BF'` cf. Equation (8) [[1]](#1), or `coeff='U'` cf. Equation (14) [[1]](#1), depending on the coefficient or set of coefficients required from Table 2(a) or Table 2(b).  See docstring for different implementations of this method.
* `get_B` : Method may take integer values of `2`, `4`, or `6` as arguments for *k*; these results in Table 1 should agree with Equation (6) [[1]](#1).
* `yamazaki2file` : <*table*> should be given as a string argument and entered as `'T1'` for Table 1, `'T2A'` for Table 2(a), or `'T2B'` for Table 2(b) [[1]](#1); <*format*> should also be given as a string argument and entered as `'CSV'` or `'JSON'`.
* `rosebrink2file` : <*table*> should be given as a string argument as explained in the corresponding docstring by calling `help` on the method; <*format*> should also be given as a string argument and entered as `'CSV'` or `'JSON'`.


| Returned quantity | Method | Arguments |
| --- | --- | --- |
| $\texttt{DataFrame}$ representation of Table 1 [[1]](#1) | `get_table1` | <sup>*</sup><*m*> |
| $\texttt{DataFrame}$ representation of Table 2(a) [[1]](#1) | `get_table2a` | *None* |
| $\texttt{DataFrame}$ representation of Table 2(b) [[1]](#1) | `get_table2b` | *None* |
| $\texttt{List}$ corresponding to data in specified row of Table 1 [[1]](#1) | `get_row_table1` | *J* |
| $\texttt{List}$ or $\texttt{float}$ corresponding to data in specified row of Table 2(a) or 2(b) [[1]](#1) depending on input arguments | `get_row_table2` | *J<sub>i</sub>*, *J<sub>f</sub>*, *L<sub>1</sub>*, *L<sub>2</sub>*, <sup>*</sup>*k*, *<sup>**</sup>coeff* |
| $\texttt{Float}$ corresponding to *B<sub>k</sub>(J)* listed in Table 1 [[1]](#1) | `get_B` | *k*, *J* |
| Dumps specified Table from Yamazaki [[1]](#1) to file in current working directory in a CSV or JSON format | `yamazaki2file` | *table*, *format* |
| $\texttt{DataFrame}$ representation of $R_{k}(L_{1}L_{2}J_{i}J_{f})$ table (integral $J$) [[2]](#2) | `get_tableRa` | *None* |
| $\texttt{DataFrame}$ representation of  $R_{k}(L_{1}L_{2}J_{i}J_{f})$ table (half-integral $J$) [[2]](#2) | `get_tableRb` | *None* |
| $\texttt{DataFrame}$ representation of $U_{k}(L_{1} J_{i} J_{f})$ and $U_{k}(L_{2} J_{i} J_{f})$ table $(L_{2} = L_{1} + 1)$ [[2]](#2) | `get_tableU` | *None* |
| $\texttt{DataFrame}$ representation of $S_{k}(l_{1}l_{2}J s)$ table (integral $J$) [[2]](#2) | `get_tableSa` | *None* |
| $\texttt{DataFrame}$ representation of $S_{k}(l_{1}l_{2}J s)$ table (half-integral $J$) [[2]](#2) | `get_tableSb` | *None* |
| $\texttt{DataFrame}$ representation of $\rho_{k}(J m)$ table (integral $J$) [[2]](#2) | `get_tablePa` | *None* |
| $\texttt{DataFrame}$ representation of $\rho_{k}(J m)$ table (half-integral $J$) [[2]](#2) | `get_tablePb` | *None* |
| Dumps specified Table from Rose and Brink [[2]](#2) to file in current working directory in a CSV or JSON format | `rosebrink2file` | *table*, *format* |


## References

<a id="1">[1]</a>
T. Yamazaki,
*"Tables of Coefficients for Angular Distribution of Gamma Rays from Aligned Nuclei"*,
Nucl. Data, Sect. A, Vol. **3**, Num. 1 (1967);
https://doi.org/10.1016/S0550-306X(67)80002-8.

<a id="2">[2]</a>
H.J. Rose, D.M. Brink,
*"Angular Distributions of Gamma Rays in Terms of Phase-Defined Reduced Matrix Elements"*,
Rev. Mod. Phys., Vol. **39**, Num. 2, p. 306 (1967);
https://doi.org/10.1103/RevModPhys.39.306.