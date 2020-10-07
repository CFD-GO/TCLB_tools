# TCLB tools

[![CircleCI](https://circleci.com/gh/ggruszczynski/TCLB_tools.svg?style=shield&circle-token=:circle-token)](https://circleci.com/gh/ggruszczynski/TCLB_tools)
[![codecov](https://codecov.io/gh/ggruszczynski/TCLB_tools/badge.svg)](https://codecov.io/gh/ggruszczynski/TCLB_tools)


Various tools for the [TCLB](/CFD-GO/TCLB/) project.

## Structure
All tools are divided by language: [R](R/), [Python](Python/), etc. In each directory, the subdirectories should be specific packages related to different aspects of interaction with TCLB

## Installation
### R

Supervisor: [Łukasz Łaniewski-Wołłk](https://github.com/llaniewski)

```R
devtools::install_github("CFD-GO/TCLB_tools/R/TCLBtools")
```

### Python

Supervisor: [Michał Dzikowski](https://github.com/mdzik)

Content of Python/ directory is meant as a packages repository. Add XXXX/Python to your PYTHONPATH than import packages by subdirectory name.

### Python: symbolic_tools

Supervisor: [Grzegorz Gruszczyński](https://github.com/ggruszczynski)

It is recomended to use PyCharm to run the scripts.
It facilitates recognition of paths and folders.
Create a PyCharm project by opening the `symbolic_tools/` directory in the editor.
Run the `examples/` . Some are listed below.

#### Example: Forcing terms

```.py
from SymbolicCollisions.core.printers import print_as_vector
from SymbolicCollisions.core.ContinousCMTransforms import \
    ContinousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D
from SymbolicCollisions.core.cm_symbols import rho, moments_dict


ccmt = ContinousCMTransforms(dzeta3D, u3D, F3D, rho)
lattice = 'D2Q9'

print("\n--- FORCES ---")
print('\n\n// === continuous central moments === \n ')

print('\n//Force -> Force_cm - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
F_cm = get_mom_vector_from_continuous_def(ccmt.get_force_He_MB,
                                          continuous_transformation=ccmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(F_cm, 'F_cm')

```

Output:
```.sh
//Force -> Force_cm - from continous definition: 
k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) 
where fun = forceM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o 
	F_cm[0] = 0;
	F_cm[1] = Fhydro.x*m00/rho;
	F_cm[2] = Fhydro.y*m00/rho;
	F_cm[3] = 0;
	F_cm[4] = 0;
	F_cm[5] = 0;
	F_cm[6] = 1/3.*Fhydro.y*m00/rho;
	F_cm[7] = 1/3.*Fhydro.x*m00/rho;
	F_cm[8] = 0;
```

#### Example: Equilibrium Distribution

```.py
from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix
from sympy import Symbol
from SymbolicCollisions.core.ContinousCMTransforms import ContinousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D

from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict

from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def


lattice = 'D2Q9'
ccmt = ContinousCMTransforms(dzeta3D, u3D, F3D, rho)
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)

print("\n--- EQUILIBRIA ---")
print('// === discrete cm ===\n ')

print('\n//population_eq -> cm_eq - by definition: k_mn = sum( (e_ix-ux)^m (e_iy-uy)^n * population_eq_i)')
print("moments: first order (linear) velocity expansion.")
pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma_first_order(i),
                                          discrete_transform=dcmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(pop_eq, 'cm_eq_first_order')

print("moments: second order (quadratic) velocity expansion.")
pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma(i),
                                          discrete_transform=dcmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(pop_eq, 'cm_eq_second_order')

print('\n\n// === continous cm === \n ')
# to calculate particular moment
row = moments_dict['D2Q9'][0]
moment = ccmt.get_cm(row, ccmt.get_Maxwellian_DF)
print_as_vector(Matrix([moment]), 'particular_moment')


print('\n//population_eq -> cm_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])
print_as_vector(cm_eq, 'cm_eq')
```

Output
```.sh
--- EQUILIBRIA ---
// === discrete cm ===
 

//population_eq -> cm_eq - by definition: k_mn = sum( (e_ix-ux)^m (e_iy-uy)^n * population_eq_i)
moments: first order (linear) velocity expansion.
	cm_eq_first_order[0] = m00;
	cm_eq_first_order[1] = 0;
	cm_eq_first_order[2] = 0;
	cm_eq_first_order[3] = m00*(-ux2 + 1/3.);
	cm_eq_first_order[4] = m00*(-uy2 + 1/3.);
	cm_eq_first_order[5] = -m00*uxuy;
	cm_eq_first_order[6] = 2.*m00*ux2*u.y;
	cm_eq_first_order[7] = 2.*m00*u.x*uy2;
	cm_eq_first_order[8] = m00*(-3.*ux2*uy2 - 1/3.*ux2 - 1/3.*uy2 + 1/9.);
moments: second order (quadratic) velocity expansion.
	cm_eq_second_order[0] = m00;
	cm_eq_second_order[1] = 0;
	cm_eq_second_order[2] = 0;
	cm_eq_second_order[3] = 1/3.*m00;
	cm_eq_second_order[4] = 1/3.*m00;
	cm_eq_second_order[5] = 0;
	cm_eq_second_order[6] = -m00*ux2*u.y;
	cm_eq_second_order[7] = -m00*u.x*uy2;
	cm_eq_second_order[8] = m00*(3.*ux2*uy2 + 1/9.);


// === continous cm === 
 
	particular_moment[0] = m00;

//population_eq -> cm_eq - from continous definition: 
k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) 
where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o 
	cm_eq[0] = m00;
	cm_eq[1] = 0;
	cm_eq[2] = 0;
	cm_eq[3] = 1/3.*m00;
	cm_eq[4] = 1/3.*m00;
	cm_eq[5] = 0;
	cm_eq[6] = 0;
	cm_eq[7] = 0;
	cm_eq[8] = 1/9.*m00;

```

#### Example: Collision Kernel

```.py
from sympy.matrices import eye
from sympy.printing import print_ccode
from SymbolicCollisions.core.cm_symbols import omega_ade, omega_b, omega_v, m00
from SymbolicCollisions.core.cm_symbols import dynamic_import
from SymbolicCollisions.core.DiscreteCMTransforms import get_DF, get_m00
from SymbolicCollisions.core.printers import print_u2, print_as_vector
from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix, get_shift_matrix

# inspired by:
# "Consistent Forcing Scheme in the cascaded LBM" L. Fei et al. 2017
# eqs 8-12 : (eye(q)-S)*cm + S*cm_eq + (eye(q)-S/2.)*force_in_cm_space

# SETUP
d = 2
q = 9
model = 'hydro'  # choose from '['hydro', 'ade', 'ade_with_f']

# DYNAMIC IMPORTS
ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")
if d == 3:
    ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
else:
    ez = None


def get_s_relax_switcher(choice):
    s_relax_switcher = {
        'hydro': ("SymbolicCollisions.core.cm_symbols", f"S_relax_hydro_D{d}Q{q}"),
        'ade_with_f': ("SymbolicCollisions.core.cm_symbols", f"S_relax_ADE_D{d}Q{q}"),
        'ade': ("SymbolicCollisions.core.cm_symbols", f"S_relax_ADE_D{d}Q{q}"),
    }
    which_model = s_relax_switcher.get(choice, lambda: "Invalid argument")
    return dynamic_import(*which_model)


S_Relax = get_s_relax_switcher(model)

hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}")
hardcoded_F_cm = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_pf_D{d}Q{q}")
from SymbolicCollisions.core.cm_symbols import Force_str as F_str

# ARRANGE STUFF
Mraw = get_raw_moments_matrix(ex, ey, ez)
Nraw = get_shift_matrix(Mraw.inv(), ex, ey, ez)

# from sympy import pprint
# pprint(Mraw)  # see what you have done
# pprint(Nraw)

pop_in_str = 'x_in'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations
cm_eq_pop_str = 'cm_eq'  # symbol defining populations


# GENERATE CODE
def make_header(choice):
    model_switcher = {
        'hydro': f"CudaDeviceFunction void relax_and_collide_hydro_with_F(real_t {pop_in_str}[{q}], real_t {omega_v}, vector_t u, vector_t {F_str}) \n{{",
        'ade_with_f': f"CudaDeviceFunction void relax_and_collide_ADE_with_F(real_t {pop_in_str}[{q}], real_t {omega_ade}, vector_t u, vector_t {F_str}) \n{{",
        'ade': f"CudaDeviceFunction void relax_and_collide_ADE(real_t {pop_in_str}[{q}], real_t {omega_ade}, vector_t u) \n{{",
    }
    result = model_switcher.get(choice, lambda: "Invalid argument")
    print(result)


make_header(model)

print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
# print(f"real_t {sv} = omega;")
# print("real_t bulk_visc = 1./6. ;")
# print("real_t {sb} = 1./(3*bulk_visc + 0.5);")
# print(f"real_t {sb} = omega_bulk;\n")  # s_b = 1./(3*bulk_visc + 0.5)
print_u2(d)
print_ccode(get_m00(q, pop_in_str), assign_to=f'\treal_t {m00}')


def make_variables(choice):
    model_switcher = {
        'hydro': f"\n\treal_t {temp_pop_str}[{q}];\n",
        'ade_with_f': f"\n\treal_t {temp_pop_str}[{q}];\n",
        'ade': f"\n\treal_t {temp_pop_str}[{q}];\n",
    }
    # Get the function from switcher dictionary
    result = model_switcher.get(choice, lambda: "Invalid argument")
    print(result)


make_variables(model)

print(f"\tfor (int i = 0; i < {q}; i++) {{\n\t"
      f"\t{temp_pop_str}[i] = {pop_in_str}[i];}}")

populations = get_DF(q, pop_in_str)
temp_populations = get_DF(q, temp_pop_str)
cm_eq = get_DF(q, cm_eq_pop_str)
F_cm = get_DF(q, F_str)
m = Mraw * temp_populations

print("\n\t//raw moments from density-probability functions")
# print("\t//[m00, m10, m01, m20, m02, m11, m21, m12, m22]")
print_as_vector(m, print_symbol=pop_in_str)

print("\n\t//central moments from raw moments")
cm = Nraw * populations
print_as_vector(cm, print_symbol=temp_pop_str)

print("\n\t//collision in central moments space")
# print("//calculate equilibrium distributions in cm space")
# print("real_t {cm_eq_pop_str}[{q}];\n")
# print_as_vector(hardcoded_cm_eq, cm_eq_pop_str)  # save time, verbosity
# print("//calculate forces in cm space")
# print("real_t {F_cm_str}[{q}];")
# print_as_vector(hardcoded_F_cm, F_cm_str)  # save time, verbosity
print("\t//collide")


def make_collision(choice):
    model_switcher = {
        # Relax 2nd moments for hydro, SOI
        'hydro': (eye(q) - S_Relax) * temp_populations
                 + S_Relax * hardcoded_cm_eq
                 + (eye(q) - S_Relax / 2) * hardcoded_F_cm,

        # Relax 1st moments for ADE, SOI
        'ade_with_f': (eye(q) - S_Relax) * temp_populations
                      + S_Relax * hardcoded_cm_eq
                      + (eye(q) - S_Relax / 2) * hardcoded_F_cm,
        # Relax 1st moments for ADE, SOI without force
        'ade': (eye(q) - S_Relax) * temp_populations
               + S_Relax * hardcoded_cm_eq,
    }
    # Get the function from switcher dictionary
    cm_after_collision = model_switcher.get(choice, lambda: "Invalid argument")
    print_as_vector(cm_after_collision, print_symbol=pop_in_str)


make_collision(model)
print("\n\t//back to raw moments")
m = Nraw.inv() * populations
print_as_vector(m, print_symbol=temp_pop_str)

print("\n\t//back to density-probability functions")
populations = Mraw.inv() * temp_populations
print_as_vector(populations, print_symbol=pop_in_str)

print("\n}\n")
```

Output
```.sh

CudaDeviceFunction void relax_and_collide_hydro_with_F(real_t x_in[9], real_t omega_nu, vector_t u, vector_t F) 
{
	//=== THIS IS AUTOMATICALLY GENERATED CODE ===
	real_t uxuy = u.x*u.y;
	real_t ux2 = u.x*u.x;
	real_t uy2 = u.y*u.y;

real_t m00 = x_in[0] + x_in[1] + x_in[2] + x_in[3] + x_in[4] + x_in[5] + x_in[6] + x_in[7] + x_in[8];

	real_t temp[9];

	for (int i = 0; i < 9; i++) {
		temp[i] = x_in[i];}

	//raw moments from density-probability functions
	x_in[0] = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7] + temp[8];
	x_in[1] = temp[1] - temp[3] + temp[5] - temp[6] - temp[7] + temp[8];
	x_in[2] = temp[2] - temp[4] + temp[5] + temp[6] - temp[7] - temp[8];
	x_in[3] = temp[1] + temp[3] + temp[5] + temp[6] + temp[7] + temp[8];
	x_in[4] = temp[2] + temp[4] + temp[5] + temp[6] + temp[7] + temp[8];
	x_in[5] = temp[5] - temp[6] + temp[7] - temp[8];
	x_in[6] = temp[5] + temp[6] - temp[7] - temp[8];
	x_in[7] = temp[5] - temp[6] - temp[7] + temp[8];
	x_in[8] = temp[5] + temp[6] + temp[7] + temp[8];

	//central moments from raw moments
	temp[0] = x_in[0];
	temp[1] = -u.x*x_in[0] + x_in[1];
	temp[2] = -u.y*x_in[0] + x_in[2];
	temp[3] = ux2*x_in[0] - 2.*u.x*x_in[1] + x_in[3];
	temp[4] = uy2*x_in[0] - 2.*u.y*x_in[2] + x_in[4];
	temp[5] = uxuy*x_in[0] - u.x*x_in[2] - u.y*x_in[1] + x_in[5];
	temp[6] = -ux2*u.y*x_in[0] + ux2*x_in[2] + 2.*uxuy*x_in[1] - 2.*u.x*x_in[5] - u.y*x_in[3] + x_in[6];
	temp[7] = -u.x*uy2*x_in[0] + 2.*uxuy*x_in[2] - u.x*x_in[4] + uy2*x_in[1] - 2.*u.y*x_in[5] + x_in[7];
	temp[8] = ux2*uy2*x_in[0] - 2.*ux2*u.y*x_in[2] + ux2*x_in[4] - 2.*u.x*uy2*x_in[1] + 4.*uxuy*x_in[5] - 2.*u.x*x_in[7] + uy2*x_in[3] - 2.*u.y*x_in[6] + x_in[8];

	//collision in central moments space
	//collide
	x_in[0] = m00;
	x_in[1] = 1/2.*F.x;
	x_in[2] = 1/2.*F.y;
	x_in[3] = 1/3.*m00*omega_bulk - 1/2.*omega_bulk*temp[3] - 1/2.*omega_bulk*temp[4] - 1/2.*omega_nu*temp[3] + 1/2.*omega_nu*temp[4] + temp[3];
	x_in[4] = 1/3.*m00*omega_bulk - 1/2.*omega_bulk*temp[3] - 1/2.*omega_bulk*temp[4] + 1/2.*omega_nu*temp[3] - 1/2.*omega_nu*temp[4] + temp[4];
	x_in[5] = -temp[5]*(omega_nu - 1.);
	x_in[6] = 1/6.*F.y;
	x_in[7] = 1/6.*F.x;
	x_in[8] = 1/9.*m00;

	//back to raw moments
	temp[0] = x_in[0];
	temp[1] = u.x*x_in[0] + x_in[1];
	temp[2] = u.y*x_in[0] + x_in[2];
	temp[3] = ux2*x_in[0] + 2.*u.x*x_in[1] + x_in[3];
	temp[4] = uy2*x_in[0] + 2.*u.y*x_in[2] + x_in[4];
	temp[5] = uxuy*x_in[0] + u.x*x_in[2] + u.y*x_in[1] + x_in[5];
	temp[6] = ux2*u.y*x_in[0] + ux2*x_in[2] + 2.*uxuy*x_in[1] + 2.*u.x*x_in[5] + u.y*x_in[3] + x_in[6];
	temp[7] = u.x*uy2*x_in[0] + 2.*uxuy*x_in[2] + u.x*x_in[4] + uy2*x_in[1] + 2.*u.y*x_in[5] + x_in[7];
	temp[8] = ux2*uy2*x_in[0] + 2.*ux2*u.y*x_in[2] + ux2*x_in[4] + 2.*u.x*uy2*x_in[1] + 4.*uxuy*x_in[5] + 2.*u.x*x_in[7] + uy2*x_in[3] + 2.*u.y*x_in[6] + x_in[8];

	//back to density-probability functions
	x_in[0] = temp[0] - temp[3] - temp[4] + temp[8];
	x_in[1] = 1/2.*temp[1] + 1/2.*temp[3] - 1/2.*temp[7] - 1/2.*temp[8];
	x_in[2] = 1/2.*temp[2] + 1/2.*temp[4] - 1/2.*temp[6] - 1/2.*temp[8];
	x_in[3] = -1/2.*temp[1] + 1/2.*temp[3] + 1/2.*temp[7] - 1/2.*temp[8];
	x_in[4] = -1/2.*temp[2] + 1/2.*temp[4] + 1/2.*temp[6] - 1/2.*temp[8];
	x_in[5] = 1/4.*temp[5] + 1/4.*temp[6] + 1/4.*temp[7] + 1/4.*temp[8];
	x_in[6] = -1/4.*temp[5] + 1/4.*temp[6] - 1/4.*temp[7] + 1/4.*temp[8];
	x_in[7] = 1/4.*temp[5] - 1/4.*temp[6] - 1/4.*temp[7] + 1/4.*temp[8];
	x_in[8] = -1/4.*temp[5] - 1/4.*temp[6] + 1/4.*temp[7] + 1/4.*temp[8];

}
```