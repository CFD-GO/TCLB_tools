**Paper under review**

# On recovering the second-order convergence of the lattice Boltzmann method with reaction-type source terms

July 2021
Grzegorz Gruszczyński, Michał Dzikowski, Łukasz Łaniewski-Wołłk

<https://arxiv.org/abs/2107.03962>

- run* sciripts are intended to test varius variants of convergence of the investigated schemes
- AC is Allen-Cahn abbreviation
- TCLB model is in models/reaction/d2q9_AllenCahn_SourceTerm_SOI

## 3.3. The Allen-Cahn equation - Solidification example

Responsible for benchmark scripts
@ggruszczynski
@mdzik

The source term solver script is here:
`Python/TCLB_tools/SymbolicCollisions/examples/lb_code_generic_generators/allen_cahn_source_term_soi/source_term_solver.py`

```.sh
$ ./configure --enable-rinside --enable-double --enable-graphics --enable-cpp11 --with-cuda-arch=sm_75
$ make -j8 d2q9_AllenCahn_SourceTerm_SOI

# run examples
$ CLB/d2q9_AllenCahn_SourceTerm_SOI/main example/reaction/d2q9_AllenCahn_SourceTerm_SOI.xml 
$ CLB/d2q9_AllenCahn_SourceTerm_SOI example/experimental/d2q9_AllenCahn_SourceTerm_SOI_runr.xml
```

## 3.3.3 Periodic 2D benchmark — impact of the Damköhler number on convergence

Responsible for benchmark scripts
@ggruszczynski

```.sh
# run the 'reverse' convergence study: from dense to sparse lattice
$ python3 run2D_AC_scalling_dense2sparse.py # to run cases, save results as dataframe.pkl
$ python3 proces_df_dense2sparse # to process the results

# alternatively:
# run the convergence study: from sparse to dense lattice
$ python3 run2D_AC_scalling_sparse2dense.py
$ python3 process_df_sparse2dense
```
