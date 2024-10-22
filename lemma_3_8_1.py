# This script is designed to verify the correctness of the calculations in Lemma 3.8 (Case III (i)) of the paper,
# allowing reviewers and readers to easily check and validate the results.
# The calculations are complex, involving judgments of polynomial parity and signs across multiple parameters.
# Including the full details in the paper would take up significant space and make verification more difficult.

import sympy as sp

# Define symbols
x = sp.symbols('x')
m, t, k, a = sp.symbols('m t k a', real=True)
m_lower_bound = 3*t + 2*k + 3
t_lower_bound = 3
k_lower_bound = 2
a_lower_bound = 0

# Define matrix: matrix_lem_3_8_1 (M_3)
matrix_lem_3_8_1 = sp.Matrix([
    [1, t, 1, k],
    [2, 0, 1, 0],
    [2, t, 0, 0],
    [2, 0, 0, 0]
])
# The characteristic polynomial of matrix_lem_3_8_1
char_poly = matrix_lem_3_8_1.charpoly(x).as_expr()

# Define the expression for x
value_to_substitute = sp.sqrt((1+sp.sqrt(m-3))**2 + a - (m - 3*t - 2*k - 3))

# The formular of \Phi_{M_3}\left( \sqrt{(1+\sqrt{m-3})^2+a-(m-3t-2k-3)}\right)
poly_substituted = char_poly.subs(x, value_to_substitute)

# Define psi_1 and psi_2 expressions
psi_1 = (4*a + 4*k + 6*t)*sp.sqrt(m - 3) + 4*m - 2*k - 3*t + 2*a*k + 3*a*t + 2*k*t + a**2 - 13

psi_2 = (a + 2*k + 6*t + 2*sp.sqrt(m - 3) + 1) * sp.sqrt(a + 2*k + 3*t + 2*sp.sqrt(m - 3) + 1)

#Compare
result_1 = sp.simplify(poly_substituted - (psi_1 - psi_2))

# Output
print(f'The characteristic polynomial of matrix_lem_3_8_1: {char_poly}')
if result_1 == 0:
    print("\nIs the formula of Phi_{M_3}(value_to_substitute) in paper correct?: correct")
else:
    print("\nIs the formula of Phi_{M_3}(value_to_substitute) in paper correct?: error")

"""The following verifies that the calculations in the paper are correct."""
# Compute psi_1^2(m,l,r,a) - psi_2^2(m,l,r,a)
diff_psi_squared = sp.expand(psi_1**2 - psi_2**2)
#The formula of psi_1^2 - psi_2^2 given in paper
diff_psi_squared_coefficients_1 = (
     (32*a + 32*k + 48*t - 8)*m + 8*a**3 + 24*a**2*k + 36*a**2*t - 6*a**2
     + 16*a*k**2 + 64*a*k*t - 40*a*k + 36*a*t**2 - 84*a*t - 116*a + 16*k**2*t - 40*k**2
     + 24*k*t**2 - 168*k*t - 128*k - 180*t**2 - 216*t + 18
)

diff_psi_squared_coefficients_2 = (
    16*m**2 + (24*a**2 + 48*a*k + 72*a*t - 12*a + 16*k**2 + 64*k*t - 40*k + 36*t**2
    - 84*t - 116)*m + a**4 + 4*a**3*k + 6*a**3*t - a**3 + 4*a**2*k**2 + 16*a**2*k*t - 10*a**2*k
    + 9*a**2*t**2 - 21*a**2*t - 77*a**2 + 8*a*k**2*t - 20*a*k**2 + 12*a*k*t**2 - 84*a*k*t
    - 160*a*k - 90*a*t**2 - 252*a*t + 33*a - 8*k**3 + 4*k**2*t**2 - 68*k**2*t - 56*k**2 - 156*k*t**2
    - 244*k*t + 118*k - 108*t**3 - 171*t**2 + 243*t + 204
)

diff_psi_squared_in_paper = (diff_psi_squared_coefficients_1*sp.sqrt(m - 3) + diff_psi_squared_coefficients_2)
result_2 = sp.simplify(diff_psi_squared - diff_psi_squared_in_paper)
if result_2 == 0:
    print("\nIs the formula of psi_1^2 - psi_2^2 in paper correct?: correct")
else:
    print("\nIs the formula of psi_1^2 - psi_2^2 in paper correct?: error")
#The above result will show that diff_psi_squared = diff_psi_squared_in_paper
# - Both diff_psi_squared_coefficients_1 and diff_psi_squared_coefficients_2 can be regarded as functions of m.
# - We need to prove that both are monotonically increasing with respect to m.
# - diff_psi_squared_coefficients_1 is a linear function of m
# - diff_psi_squared_coefficients_2 is a quadratic function of m.
# - Taking the first derivative of diff_psi_squared_coefficients_2 with respect to m, the derivative need to be greater than 0.

#Easily, we know that the diff_psi_squared_coefficients_1 is monotonically increasing with respect to m when m > 3*t + 2*k + 3.
#It is easy to see that the value of the following polynomial is greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_1 = sp.collect(sp.simplify(diff_psi_squared_coefficients_1.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_1: {substitute_lb_into_diff_psi_squared_coefficients_1}')

diff_psi_squared_coefficients_2_derivative_m = sp.diff(diff_psi_squared_coefficients_2, m)
#From the following equation, we know that the diff_psi_squared_coefficients_2 is monotonically increasing with respect to m when m >= 3*t + 2*k + 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m = sp.collect(sp.simplify(diff_psi_squared_coefficients_2_derivative_m.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m: {substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m}')

#----------------------------------------------------------------------------------------------------------------
#Substitute m_lower_bound = 3*t + 2*k + 3 into psi_1^2 - psi_2^2.
substitute_m_into_diff_psi_squared = diff_psi_squared.subs(m,m_lower_bound)
#The formula of substitute_m_into_diff_psi_squared given in paper
substitute_m_into_coefficients_1 = ((36*a + 24*k - 36)*t**2
     + (36*a**2 + 64*a*k + 12*a + 16*k**2 + 24*k - 96)*t
     + 8*a**3 + 24*a**2*k - 6*a**2 + 16*a*k**2 + 24*a*k - 20*a
     + 24*k**2 - 48*k - 6
)

substitute_m_into_coefficients_2 = (
    ((9*a**2 + 12*a*k + 126*a + 4*k**2 + 108*k - 171)*t**2 + (6*a**3 + 16*a**2*k
     + 51*a**2 + 8*a*k**2 + 204*a*k - 72*a + 108*k**2 - 148*k - 69)*t + a**4 + 4*a**3*k - a**3
     + 4*a**2*k**2 + 38*a**2*k - 5*a**2) + 76*a*k**2 - 40*a*k - 3*a + 24*k**3 - 24*k**2 - 42*k
)


substitute_m_into_diff_psi_squared_in_paper = (substitute_m_into_coefficients_1*sp.sqrt(3*t+2*k) + substitute_m_into_coefficients_2)
# The above result will show that substitute_m_into_diff_psi_squared = substitute_m_into_diff_psi_squared_in_paper
# Both substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 can be regarded as functions of t.
# We need to show that both are monotonically increasing with respect to t.
# It is evident that substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2
# are both quadratic functions of t.
# If the first derivatives of substitute_m_into_coefficients_1
# and substitute_m_into_coefficients_2 with respect to t are all greater than 0,
# then substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 will be monotonically increasing with respect to t.
result_3 = sp.simplify(substitute_m_into_diff_psi_squared - substitute_m_into_diff_psi_squared_in_paper)
if result_3 == 0:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: error")

substitute_m_into_coefficients_1_derivative_t_1 = sp.diff(substitute_m_into_coefficients_1, t)
#From the following equation, we know that the substitute_m_into_coefficients_1 is monotonically increasing with respect to t when t >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_1_derivative_t_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_t_1.subs(t, t_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_t_1: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_t_1}')

substitute_m_into_coefficients_2_derivative_t_1 = sp.diff(substitute_m_into_coefficients_2, t)
#From the following equation, we know that the substitute_m_into_coefficients_2 is monotonically increasing with respect to t when t >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_2_derivative_t_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_t_1.subs(t, t_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_t_1: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_t_1}')

#--------------------------------------------------------------------------------------------------------------------
#Substitute t_lower_bound = 3 into substitute_m_into_diff_psi_squared.
substitute_m_t_into_diff_psi_squared = substitute_m_into_diff_psi_squared.subs(t,t_lower_bound)
#The formula of substitute_m_t_into_diff_psi_squared given in paper
substitute_m_t_into_coefficients_1 = (
     (16*a + 72)*k**2 + (24*a**2 + 216*a + 240)*k + 8*a**3 + 102*a**2 + 340*a - 618
)

substitute_m_t_into_coefficients_2 = (
    24*k**3 + (4*a**2 + 100*a + 336)*k**2 + (4*a**3 + 86*a**2 + 680*a + 486)*k
    + a**4 + 17*a**3 + 229*a**2 + 915*a - 1746
)

substitute_m_t_into_diff_psi_squared_in_paper = (substitute_m_t_into_coefficients_1*sp.sqrt(2*k + 9)+substitute_m_t_into_coefficients_2)
result_4 = sp.simplify(substitute_m_t_into_diff_psi_squared - substitute_m_t_into_diff_psi_squared_in_paper)
if result_4 == 0:
    print("\nIs the formula of substitute_m_t_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_t_into_diff_psi_squared in paper correct?: error")
# The above result will show that substitute_m_t_into_diff_psi_squared = substitute_m_t_into_diff_psi_squared_in_paper
# Both substitute_m_t_into_coefficients_1 and substitute_m_t_into_coefficients_2 can be regarded as functions of k.
# We need to show that both are monotonically increasing with respect to k.
# It is evident that substitute_m_t_into_coefficients_1 is a quadratic function of k,
# and substitute_m_t_into_coefficients_2 is a cubic functions of k.
# If the first derivative of substitute_m_t_into_coefficients_1 with respect to k,
# and the first and second derivatives of substitute_m_t_into_coefficients_2 are all greater than 0,
# then substitute_m_t_into_coefficients_1 and substitute_m_t_into_coefficients_2 will be monotonically increasing with respect to k.
substitute_m_t_into_coefficients_1_derivative_k_1 = sp.diff(substitute_m_t_into_coefficients_1, k)
#From the following equation, we know that the substitute_m_t_into_coefficients_1 is monotonically increasing with respect to k when k >= 2.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_t_into_coefficients_1_derivative_k_1 = sp.collect(sp.simplify(substitute_m_t_into_coefficients_1_derivative_k_1.subs(k, k_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_t_into_coefficients_1_derivative_k_1: {substitute_lb_into_substitute_m_t_into_coefficients_1_derivative_k_1}')


substitute_m_t_into_coefficients_2_derivative_k_1 = sp.diff(substitute_m_t_into_coefficients_2, k)
substitute_m_t_into_coefficients_2_derivative_k_2 = sp.diff(substitute_m_t_into_coefficients_2_derivative_k_1, k)
#From the following equation, we know that the substitute_m_t_into_coefficients_2 is monotonically increasing with respect to r when k >= 2.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_k_1 = sp.collect(sp.simplify(substitute_m_t_into_coefficients_2_derivative_k_1.subs(k, k_lower_bound)),a)
substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_k_2 = sp.collect(sp.simplify(substitute_m_t_into_coefficients_2_derivative_k_2.subs(k, k_lower_bound)),a)
print(f'\nThe first derivative function of substitute_m_t_into_coefficients_2 with respect to k: {sp.collect(sp.expand(substitute_m_t_into_coefficients_2_derivative_k_1),k)}')
print(f'The value of substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_k_1: {substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_k_1}')
print(f'The second derivative function of substitute_m_t_into_coefficients_2 with respect to k: {sp.collect(sp.expand(substitute_m_t_into_coefficients_2_derivative_k_2),k)}')
print(f'The value of substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_k_2: {substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_k_2}')

#--------------------------------------------------------------------------------------------------------------
#Substitute k_lower_bound = 2 into substitute_m_t_into_diff_psi_squared.
substitute_m_t_k_into_diff_psi_squared = substitute_m_t_into_diff_psi_squared.subs(k,k_lower_bound)
#The formula of substitute_m_t_k_into_diff_psi_squared given in paper
substitute_m_t_k_into_diff_psi_squared_in_paper =(
        2675 * a + 836 * sp.sqrt(13) * a + 150 * sp.sqrt(13)
        + 150 * sp.sqrt(13) * a ** 2 + 8 * sp.sqrt(13) * a ** 3
        + 417 * a ** 2 + 25 * a ** 3 + a ** 4 + 762
)
result_5 = sp.simplify(substitute_m_t_k_into_diff_psi_squared - substitute_m_t_k_into_diff_psi_squared_in_paper)
if result_5 == 0:
    print("\nIs the formula of substitute_m_t_k_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_t_k_into_diff_psi_squared in paper correct?: error")