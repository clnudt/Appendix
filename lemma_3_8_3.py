# This script is designed to verify the correctness of the calculations in Lemma 3.8 (Case IV) of the paper,
# allowing reviewers and readers to easily check and validate the results.
# The calculations are complex, involving judgments of polynomial parity and signs across multiple parameters.
# Including the full details in the paper would take up significant space and make verification more difficult.

import sympy as sp

# Define symbols
x = sp.symbols('x')
m, p, t, a = sp.symbols('m p t a', real=True)
m_lower_bound = 3*t + p + 3
p_lower_bound = 2
t_lower_bound = 3
a_lower_bound = 0

# Define matrix: matrix_lem_3_8_3 (M_5)
matrix_lem_3_8_3 = sp.Matrix([
    [0, t, 2, p],
    [1, 0, 2, 0],
    [1, t, 1, 0],
    [1, 0, 0, 0]
])
# The characteristic polynomial of matrix_lem_3_8_3
char_poly = matrix_lem_3_8_3.charpoly(x).as_expr()

# Define the expression for x
value_to_substitute = sp.sqrt((1 + sp.sqrt(m - 3))**2 + a - (m - 3*t - p - 3))

# The formular of \Phi_{M_{5}}\left(\sqrt{(1+\sqrt{m-3})^2+a-(m-3t-p-3)} \right)
poly_substituted = char_poly.subs(x, value_to_substitute)

# Define psi_5 and psi_6 expressions
psi_5 = (4*a + 2*p + 6*t)*sp.sqrt(m - 3) + 4*m - p - 3*t + a*p + 3*a*t + 2*p*t + a**2 - 13

psi_6 = (a + 6*t + 2*sp.sqrt(m - 3) + 1)*sp.sqrt(a + p + 3*t + 2*sp.sqrt(m - 3) + 1)

#Compare
result_1 = sp.simplify(poly_substituted - (psi_5 - psi_6))

# Output
print(f'The characteristic polynomial of matrix_lem_3_8_3: {char_poly}')
if result_1 == 0:
    print("\nIs the formula of Phi_{M_5}(value_to_substitute) in paper correct?: correct")
else:
    print("\nIs the formula of Phi_{M_5}(value_to_substitute) in paper correct?: error")

"""The following verifies that the calculations in the paper are correct."""
# Compute psi_5^2(m,p,t,a) - psi_6^2(m,p,t,a)
diff_psi_squared = sp.expand(psi_5**2 - psi_6**2)
#The formula of psi_5^2 - psi_6^2 given in paper
diff_psi_squared_coefficients_1 = (
    (32*a + 16*p + 48*t - 8)*m + 8*a**3 + 12*a**2*p + 36*a**2*t - 6*a**2 +
    4*a*p**2 + 40*a*p*t - 12*a*p + 36*a*t**2 - 84*a*t - 116*a +
    8*p**2*t - 4*p**2 + 24*p*t**2 - 48*p*t - 56*p - 180*t**2 - 216*t + 18
)

diff_psi_squared_coefficients_2 = (
    16*m**2 + (24*a**2 + 24*a*p + 72*a*t - 12*a + 4*p**2 + 40*p*t - 12*p + 36*t**2 - 84*t - 116)*m +
    a**4 + 2*a**3*p + 6*a**3*t - a**3 + a**2*p**2 + 10*a**2*p*t - 3*a**2*p +
    9*a**2*t**2 - 21*a**2*t - 77*a**2 + 4*a*p**2*t - 2*a*p**2 + 12*a*p*t**2 - 24*a*p*t -
    76*a*p - 90*a*t**2 - 252*a*t + 33*a + 4*p**2*t**2 - 4*p**2*t - 11*p**2 -
    48*p*t**2 - 130*p*t + 37*p - 108*t**3 - 171*t**2 + 243*t + 204
)

diff_psi_squared_in_paper = (diff_psi_squared_coefficients_1*sp.sqrt(m - 3) + diff_psi_squared_coefficients_2)
result_2 = sp.simplify(diff_psi_squared - diff_psi_squared_in_paper)
if result_2 == 0:
    print("\nIs the formula of psi_5^2 - psi_6^2 in paper correct?: correct")
else:
    print("\nIs the formula of psi_5^2 - psi_6^2 in paper correct?: error")
#The above result will show that diff_psi_squared = diff_psi_squared_in_paper
# - Both diff_psi_squared_coefficients_1 and diff_psi_squared_coefficients_2 can be regarded as functions of m.
# - We need to prove that both are monotonically increasing with respect to m.
# - diff_psi_squared_coefficients_1 is a linear function of m
# - diff_psi_squared_coefficients_2 is a quadratic function of m.
# - Taking the first derivative of diff_psi_squared_coefficients_2 with respect to m, the derivative need to be greater than 0.

#Easily, we know that the diff_psi_squared_coefficients_1 is monotonically increasing with respect to m when m >= 3*t + p + 3.
#It is easy to see that the value of the following polynomial is greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_1 = sp.collect(sp.simplify(diff_psi_squared_coefficients_1.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_1: {substitute_lb_into_diff_psi_squared_coefficients_1}')

diff_psi_squared_coefficients_2_derivative_m = sp.diff(diff_psi_squared_coefficients_2, m)
#From the following equation, we know that the diff_psi_squared_coefficients_2 is monotonically increasing with respect to m when m >= 3*t + p + 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m = sp.collect(sp.simplify(diff_psi_squared_coefficients_2_derivative_m.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m: {substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m}')

#----------------------------------------------------------------------------------------------------------------
#Substitute m_lower_bound = 3*t + p + 3 into psi_5^2 - psi_6^2.
substitute_m_into_diff_psi_squared = diff_psi_squared.subs(m,m_lower_bound)
#The formula of substitute_m_into_diff_psi_squared given in paper
substitute_m_into_coefficients_1 =  (
    (4*a + 8*t + 12)*p**2 + (12*a**2 + 40*a*t + 20*a + 24*t**2 + 48*t - 16)*p +
    8*a**3 + 36*a**2*t - 6*a**2 + 36*a*t**2 + 12*a*t - 20*a - 36*t**2 - 96*t - 6
)

substitute_m_into_coefficients_2  = (
    4*p**3 +  (a**2 + 4*a*t + 22*a + 4*t**2 + 48*t + 5)*p**2 +
    (2*a**3 + 10*a**2*t + 21*a**2 + 12*a*t**2 + 120*a*t - 16*a +
    108*t**2 - 34*t - 19)*p + a**4 + 6*a**3*t - a**3 + 9*a**2*t**2 + 51*a**2*t - 5*a**2 +
    126*a*t**2 - 72*a*t - 3*a - 171*t**2 - 69*t
)

substitute_m_into_diff_psi_squared_in_paper = (substitute_m_into_coefficients_1*sp.sqrt(3*t + p) +substitute_m_into_coefficients_2)
# The above result will show that substitute_m_into_diff_psi_squared = substitute_m_into_diff_psi_squared_in_paper
# Both substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 can be regarded as functions of p.
# We need to show that both are monotonically increasing with respect to p.
# It is evident that substitute_m_into_coefficients_1 is a quadratic function of p,
# and substitute_m_into_coefficients_2 is a cubic functions of p.
# If the first derivative of substitute_m_into_coefficients_1 with respect to p,
# and the first and second derivatives of substitute_m_into_coefficients_2 are all greater than 0,
# then substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 will be monotonically increasing with respect to p.
result_3 = sp.simplify(substitute_m_into_diff_psi_squared - substitute_m_into_diff_psi_squared_in_paper)
if result_3 == 0:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: error")

substitute_m_into_coefficients_1_derivative_p_1 = sp.diff(substitute_m_into_coefficients_1, p)
#From the following equation, we know that the substitute_m_into_coefficients_1 is monotonically increasing with respect to p when p >= 2.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_p_1.subs(p, p_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_1: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_1}')

substitute_m_into_coefficients_2_derivative_p_1 = sp.diff(substitute_m_into_coefficients_2, p)
substitute_m_into_coefficients_2_derivative_p_2 = sp.diff(substitute_m_into_coefficients_2_derivative_p_1, p)
#From the following equation, we know that the substitute_m_into_coefficients_2 is monotonically increasing with respect to p when p >= 2.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_p_1.subs(p, p_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_2 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_p_2.subs(p, p_lower_bound)),a)
print(f'\nThe first derivative function of substitute_m_into_coefficients_2 with respect to p: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_p_1),p)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_1: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_1}')
print(f'The second derivative function of substitute_m_into_coefficients_2 with respect to p: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_p_2),p)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_2: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_2}')

#--------------------------------------------------------------------------------------------------------------------
#Substitute p_lower_bound = 2 into substitute_m_into_diff_psi_squared.
substitute_m_p_into_diff_psi_squared = substitute_m_into_diff_psi_squared.subs(p,p_lower_bound)
#The formula of substitute_m_p_into_diff_psi_squared given in paper
substitute_m_p_into_coefficients_1 = (
     (8*a**3 + (36*t + 18)*a**2 + (36*t**2 + 92*t + 36)*a + 12*t**2 + 32*t + 10)
)

substitute_m_p_into_coefficients_2 = (
    (9*a**2 + 150*a + 61)*t**2 + (6*a**3 + 71*a**2 + 184*a + 55)*t + a**4 + 3*a**3 + 41*a**2 + 53*a + 14
)

substitute_m_p_into_diff_psi_squared_in_paper = (substitute_m_p_into_coefficients_1*sp.sqrt(3*t + 2)+substitute_m_p_into_coefficients_2)
result_4 = sp.simplify(substitute_m_p_into_diff_psi_squared - substitute_m_p_into_diff_psi_squared_in_paper)
if result_4 == 0:
    print("\nIs the formula of substitute_m_p_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_p_into_diff_psi_squared in paper correct?: error")
# The above result will show that substitute_m_p_into_diff_psi_squared = substitute_m_p_into_diff_psi_squared_in_paper
# Both substitute_m_p_into_coefficients_1 and substitute_m_p_into_coefficients_2 can be regarded as functions of t.
# We need to show that both are monotonically increasing with respect to t.
# It is evident that substitute_m_t_into_coefficients_1 and substitute_m_t_into_coefficients_2
# are quadratic function of t.
# If the first derivative of these two polynomials with respect to t are all greater than 0,
# then substitute_m_p_into_coefficients_1 and substitute_m_p_into_coefficients_2 will be monotonically increasing with respect to t.
substitute_m_p_into_coefficients_1_derivative_t_1 = sp.diff(substitute_m_p_into_coefficients_1, t)
#From the following equation, we know that the substitute_m_p_into_coefficients_1 is monotonically increasing with respect to t when t >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_1 = sp.collect(sp.simplify(substitute_m_p_into_coefficients_1_derivative_t_1.subs(t, t_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_1: {substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_1}')


substitute_m_p_into_coefficients_2_derivative_t_1 = sp.diff(substitute_m_p_into_coefficients_2, t)
#From the following equation, we know that the substitute_m_t_into_coefficients_2 is monotonically increasing with respect to r when k >= 2.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_t_1 = sp.collect(sp.simplify(substitute_m_p_into_coefficients_2_derivative_t_1.subs(t, t_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_t_1: {substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_t_1}')

#--------------------------------------------------------------------------------------------------------------
#Substitute t_lower_bound = 3 into substitute_m_p_into_diff_psi_squared.
substitute_m_p_t_into_diff_psi_squared = substitute_m_p_into_diff_psi_squared.subs(t,t_lower_bound)
#The formula of substitute_m_p_t_into_diff_psi_squared given in paper
substitute_m_p_t_into_diff_psi_squared_in_paper =(
    1955*a + sp.sqrt(11)*(8*a**3 + 126*a**2 + 636*a + 214) +
    335*a**2 + 21*a**3 + a**4 + 728
)
result_5 = sp.simplify(substitute_m_p_t_into_diff_psi_squared - substitute_m_p_t_into_diff_psi_squared_in_paper)
if result_5 == 0:
    print("\nIs the formula of substitute_m_p_t_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_p_t_into_diff_psi_squared in paper correct?: error")