# This script is designed to verify the correctness of the calculations in Lemma 3.8 (Case III (ii)) of the paper,
# allowing reviewers and readers to easily check and validate the results.
# The calculations are complex, involving judgments of polynomial parity and signs across multiple parameters.
# Including the full details in the paper would take up significant space and make verification more difficult.

import sympy as sp

# Define symbols
x = sp.symbols('x')
m, t, p, a = sp.symbols('m t p a', real=True)
m_lower_bound = 3*t + p + 5
t_lower_bound = 3
p_lower_bound = 1
a_lower_bound = 0

# Define matrix: matrix_lem_3_8_2 (M_4)
matrix_lem_3_8_2 = sp.Matrix([
    [0, t, 1, 1, 1, p],
    [1, 0, 1, 1, 0, 0],
    [1, t, 0, 1, 1, 0],
    [1, t, 1, 0, 0, 0],
    [1, 0, 1, 0, 0, 0],
    [1, 0, 0, 0, 0, 0]
])
# The characteristic polynomial of matrix_lem_3_8_2
char_poly = matrix_lem_3_8_2.charpoly(x).as_expr()

# Define the expression for x
value_to_substitute = sp.sqrt((1 + sp.sqrt(m - 3))**2 + a - (m - 3*t - p - 5))

# The formular of \Phi_{M_{4}}\left(\sqrt{(1+\sqrt{m-3})^2+a-(m-3t-p-5)}\right)
poly_substituted = char_poly.subs(x, value_to_substitute)

# Define psi_3 and psi_4 expressions
psi_3 = (
    (6*a**2 + 8*a*p + 24*a*t + 16*a + 2*p**2 + 16*p*t + 8*p + 18*t**2 + 10*t + 8*m - 30) * sp.sqrt(m - 3)
    + (12*a + 8*p + 24*t + 16)*m
    + a**3 + 2*a**2*p + 6*a**2*t + 4*a**2 + a*p**2 + 8*a*p*t + 4*a*p + 9*a*t**2 + 5*a*t - 39*a
    + 2*p**2*t + 6*p*t**2 - 2*p*t - 30*p - 21*t**2 - 111*t - 66
)

psi_4 = (
    2*((6*t + 4)*sp.sqrt(m - 3) + 2*a + 2*p + 14*t + 3*a*t + 2*p*t + 9*t**2 + 6)
    * sp.sqrt(a + p + 3*t + 2*sp.sqrt(m - 3) + 3)
)

#Compare
result_1 = sp.simplify(poly_substituted -  (psi_3 - psi_4))

# Output
print(f'The characteristic polynomial of matrix_lem_3_8_2: {char_poly}')
if result_1 == 0:
    print("\nIs the formula of Phi_{M_4}(value_to_substitute) in paper correct?: correct")
else:
    print("\nIs the formula of Phi_{M_4}(value_to_substitute) in paper correct?: error")

"""The following verifies that the calculations in the paper are correct."""
# Compute psi_3^2(m,l,r,a) - psi_4^2(m,l,r,a)
diff_psi_squared = sp.expand(psi_3**2 - psi_4**2)
#The formula of psi_3^2 - psi_4^2 given in paper
diff_psi_squared_coefficients_1= (
    (192*a + 128*p + 384*t + 256)*m**2 +
    (160*a**3 + 320*a**2*p + 960*a**2*t + 640*a**2 + 192*a*p**2 + 1280*a*p*t + 768*a*p +
    1728*a*t**2 + 1856*a*t - 832*a + 32*p**3 + 384*p**2*t + 192*p**2 + 1152*p*t**2 + 1024*p*t -
    704*p + 864*t**3 + 432*t**2 - 3280*t - 2144)*m +
    12*a**5 + 40*a**4*p + 120*a**4*t + 80*a**4 + 48*a**3*p**2 + 320*a**3*p*t + 192*a**3*p +
    432*a**3*t**2 + 464*a**3*t - 400*a**3 + 24*a**2*p**3 + 288*a**2*p**2*t + 144*a**2*p**2 +
    864*a**2*p*t**2 + 768*a**2*p*t - 912*a**2*p + 648*a**2*t**3 + 324*a**2*t**2 - 3612*a**2*t -
    2376*a**2 + 4*a*p**4 + 96*a*p**3*t + 32*a*p**3 + 520*a*p**2*t**2 + 328*a*p**2*t - 632*a*p**2 +
    864*a*p*t**3 + 32*a*p*t**2 - 5360*a*p*t - 3072*a*p + 324*a*t**4 - 1944*a*t**3 - 10772*a*t**2 -
    10040*a*t - 348*a + 8*p**4*t + 88*p**3*t**2 + 24*p**3*t - 120*p**3 + 264*p**2*t**3 - 140*p**2*t**2 -
    1780*p**2*t - 840*p**2 + 216*p*t**4 - 1632*p*t**3 - 7832*p*t**2 - 6416*p*t + 168*p - 2700*t**4 -
    10608*t**3 - 9992*t**2 + 2364*t + 3480
)

diff_psi_squared_coefficients_2 = (
    64*m**3
    + (240*a**2 + 320*a*p + 960*a*t + 640*a + 96*p**2 + 640*p*t + 384*p + 864*t**2 + 928*t - 416)*m**2
    + (60*a**4 + 160*a**3*p + 480*a**3*t + 320*a**3 + 144*a**2*p**2 + 960*a**2*p*t + 576*a**2*p +
    1296*a**2*t**2 + 1392*a**2*t - 1200*a**2 + 48*a*p**3 + 576*a*p**2*t + 288*a*p**2 + 1728*a*p*t**2 +
    1536*a*p*t - 1824*a*p + 1296*a*t**3 + 648*a*t**2 - 7224*a*t - 4752*a + 4*p**4 + 96*p**3*t + 32*p**3 +
    520*p**2*t**2 + 328*p**2*t - 632*p**2 + 864*p*t**3 + 32*p*t**2 - 5360*p*t - 3072*p + 324*t**4 -
    1944*t**3 - 10772*t**2 - 10040*t - 348)*m
    + a**6 + 4*a**5*p + 12*a**5*t + 8*a**5 + 6*a**4*p**2 + 40*a**4*p*t + 24*a**4*p + 54*a**4*t**2 +
    58*a**4*t - 170*a**4 + 4*a**3*p**3 + 48*a**3*p**2*t + 24*a**3*p**2 + 144*a**3*p*t**2 + 128*a**3*p*t -
    472*a**3*p + 108*a**3*t**3 + 54*a**3*t**2 - 1562*a**3*t - 1036*a**3 + a**2*p**4 + 24*a**2*p**3*t +
    8*a**2*p**3 + 130*a**2*p**2*t**2 + 82*a**2*p**2*t - 446*a**2*p**2 + 216*a**2*p*t**3 + 8*a**2*p*t**2 -
    3260*a**2*p*t - 1920*a**2*p + 81*a**2*t**4 - 486*a**2*t**3 - 5285*a**2*t**2 - 5294*a**2*t + 1161*a**2 +
    4*a*p**4*t + 44*a*p**3*t**2 + 12*a*p**3*t - 156*a*p**3 + 132*a*p**2*t**3 - 70*a*p**2*t**2 -
    2042*a*p**2*t - 996*a*p**2 + 108*a*p*t**4 - 816*a*p*t**3 - 7372*a*p*t**2 - 6280*a*p*t + 2196*a*p -
    1350*a*t**4 - 7896*a*t**3 - 6292*a*t**2 + 11022*a*t + 8172*a + 4*p**4*t**2 - 12*p**4 + 24*p**3*t**3 -
    24*p**3*t**2 - 344*p**3*t - 112*p**3 + 36*p**2*t**4 - 300*p**2*t**3 - 2296*p**2*t**2 - 1496*p**2*t +
    924*p**2 - 1008*p*t**4 - 5520*p*t**3 - 3184*p*t**2 + 8940*p*t + 5544*p - 972*t**5 - 4527*t**4 +
    798*t**3 + 21153*t**2 + 20724*t + 2952
)

diff_psi_squared_in_paper = (diff_psi_squared_coefficients_1*sp.sqrt(m - 3) + diff_psi_squared_coefficients_2)
result_2 = sp.simplify(diff_psi_squared - diff_psi_squared_in_paper)
if result_2 == 0:
    print("\nIs the formula of psi_3^2 - psi_4^2 in paper correct?: correct")
else:
    print("\nIs the formula of psi_3^2 - psi_4^2 in paper correct?: error")
#The above result will show that diff_psi_squared = diff_psi_squared_in_paper
# - Both diff_psi_squared_coefficients_1 and diff_psi_squared_coefficients_2 can be regarded as functions of m.
# - We need to prove that both are monotonically increasing with respect to m.
# - It is evident that diff_psi_squared_coefficients_1 is a quadratic function of m,
# - and diff_psi_squared_coefficients_2 is a cubic function of m.
# - If the first derivative of diff_psi_squared_coefficients_1 with respect to m,
# - and the first and second derivatives of diff_psi_squared_coefficients_2 are all greater than 0,
# - then diff_psi_squared_coefficients_1 and diff_psi_squared_coefficients_2 will be monotonically increasing with respect to m.
diff_psi_squared_coefficients_1_derivative_m = sp.diff(diff_psi_squared_coefficients_1, m)
#From the following equation, we know that the diff_psi_squared_coefficients_1 is monotonically increasing with respect to m when m >= 3*t + p + 5.
#It is easy to see that the value of the following polynomial is greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_1_derivative_m_1 = sp.collect(sp.simplify(diff_psi_squared_coefficients_1_derivative_m.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_1_derivative_m_1: {substitute_lb_into_diff_psi_squared_coefficients_1_derivative_m_1}')

diff_psi_squared_coefficients_2_derivative_m_1 = sp.diff(diff_psi_squared_coefficients_2, m)
diff_psi_squared_coefficients_2_derivative_m_2 = sp.diff(diff_psi_squared_coefficients_2_derivative_m_1, m)
#From the following equation, we know that the diff_psi_squared_coefficients_2 is monotonically increasing with respect to m when m >= 3*t + p + 5.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_1 = sp.collect(sp.simplify(diff_psi_squared_coefficients_2_derivative_m_1.subs(m, m_lower_bound)),a)
substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_2 = sp.collect(sp.simplify(diff_psi_squared_coefficients_2_derivative_m_2.subs(m, m_lower_bound)),a)
print(f'\nThe first derivative function of diff_psi_squared_coefficients_2 with respect to m: {sp.collect(sp.expand(diff_psi_squared_coefficients_2_derivative_m_1),m)}')
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_1: {substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_1}')
print(f'The second derivative function of diff_psi_squared_coefficients_2 with respect to m: {sp.collect(sp.expand(diff_psi_squared_coefficients_2_derivative_m_2),m)}')
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_2: {substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_2}')

#----------------------------------------------------------------------------------------------------------------
#Substitute m_lower_bound = 3*t + 2*k + 3 into psi_3^2 - psi_4^2.
substitute_m_into_diff_psi_squared = diff_psi_squared.subs(m,m_lower_bound)
#The formula of substitute_m_into_diff_psi_squared given in paper
substitute_m_into_coefficients_1= (
    (4*a + 8*t + 32)*p**4 +
    (24*a**2 + 96*a*t + 224*a + 88*t**2 + 504*t + 360)*p**3 +
    (48*a**3 + 288*a**2*t + 464*a**2 + 520*a*t**2 + 2184*a*t + 1288*a +
    264*t**3 + 2164*t**2 + 2892*t + 952)*p**2 +
    (40*a**4 + 320*a**3*t + 352*a**3 + 864*a**2*t**2 + 2688*a**2*t +
    1328*a**2 + 864*a*t**3 + 5600*a*t**2 + 6352*a*t + 1856*a + 216*t**4 +
    2688*t**3 + 4888*t**2 + 2528*t + 264)*p +
    12*a**5 + 120*a**4*t + 80*a**4 + 432*a**3*t**2 + 944*a**3*t + 400*a**3 +
    648*a**2*t**3 + 3204*a**2*t**2 + 3108*a**2*t + 824*a**2 + 324*a*t**4 +
    3240*a*t**3 + 5164*a*t**2 + 2504*a*t + 292*a - 108*t**4 - 1536*t**3 -
    3848*t**2 - 3188*t - 840
)

substitute_m_into_coefficients_2 = (
    4*p**5 + (a**2 + 4*a*t + 48*a + 4*t**2 + 108*t + 136)*p**4
    + (4*a**3 + 24*a**2*t + 152*a**2 + 44*a*t**2 + 732*a*t + 692*a +
    24*t**3 + 784*t**2 + 1776*t + 824)*p**3
    + (6*a**4 + 48*a**3*t + 184*a**3 + 130*a**2*t**2 + 1474*a**2*t +
    1090*a**2 + 132*a*t**3 + 3386*a*t**2 + 6118*a*t + 2460*a + 36*t**4 +
    2124*t**3 + 6888*t**2 + 5976*t + 1476)*p**2
    + (4*a**5 + 40*a**4*t + 84*a**4 + 144*a**3*t**2 + 1088*a**3*t + 648*a**3 +
    216*a**2*t**3 + 4184*a**2*t**2 + 6100*a**2*t + 2160*a**2 + 108*a*t**4 +
    5664*a*t**3 + 15164*a*t**2 + 11744*a*t + 2724*a + 1908*t**4 +
    7896*t**3 + 8716*t**2 + 2948*t + 76)*p
    + a**6 + 12*a**5*t + 8*a**5 + 54*a**4*t**2 + 238*a**4*t + 130*a**4 +
    108*a**3*t**3 + 1494*a**3*t**2 + 1798*a**3*t + 564*a**3 + 81*a**2*t**4 +
    3402*a**2*t**3 + 7531*a**2*t**2 + 5266*a**2*t + 1161*a**2 +
    2538*a*t**4 + 9168*a*t**3 + 9836*a*t**2 + 3846*a*t + 412*a -
    963*t**4 - 5238*t**3 - 8491*t**2 - 5400*t - 1188
)

substitute_m_into_diff_psi_squared_in_paper = (substitute_m_into_coefficients_1*sp.sqrt(p + 3*t + 2) +substitute_m_into_coefficients_2)
# The above result will show that substitute_m_into_diff_psi_squared = substitute_m_into_diff_psi_squared_in_paper
# Both substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 can be regarded as functions of p.
# We need to show that both are monotonically increasing with respect to p.
# It is evident that substitute_m_into_coefficients_1 is a quadratic function of p
# and substitute_m_into_coefficients_2 is a quintic function of p.
# If the first, second, and third derivatives of substitute_m_into_coefficients_1 with respect to p,
# and the first, second, third and fourth derivatives of substitute_m_into_coefficients_2 with respect to p
# are all greater than 0,
# then substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 will be monotonically increasing with respect to p.
result_3 = sp.simplify(substitute_m_into_diff_psi_squared - substitute_m_into_diff_psi_squared_in_paper)
if result_3 == 0:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: error")

substitute_m_into_coefficients_1_derivative_p_1 = sp.diff(substitute_m_into_coefficients_1, p)
substitute_m_into_coefficients_1_derivative_p_2 = sp.diff(substitute_m_into_coefficients_1_derivative_p_1, p)
substitute_m_into_coefficients_1_derivative_p_3 = sp.diff(substitute_m_into_coefficients_1_derivative_p_2, p)
#From the following equation, we know that the substitute_m_into_coefficients_1 is monotonically increasing with respect to p when p >= 1.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_p_1.subs(p, p_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_2 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_p_2.subs(p, p_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_3 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_p_3.subs(p, p_lower_bound)),a)
print(f'The first derivative function of substitute_m_into_coefficients_1 with respect to p: {sp.collect(sp.expand(substitute_m_into_coefficients_1_derivative_p_1),p)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_1: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_1}')
print(f'The second derivative function of substitute_m_into_coefficients_1 with respect to p: {sp.collect(sp.expand(substitute_m_into_coefficients_1_derivative_p_2),p)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_2: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_2}')
print(f'The third derivative function of substitute_m_into_coefficients_1 with respect to p: {sp.collect(sp.expand(substitute_m_into_coefficients_1_derivative_p_3),p)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_3: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_p_3}')

substitute_m_into_coefficients_2_derivative_p_1 = sp.diff(substitute_m_into_coefficients_2, p)
substitute_m_into_coefficients_2_derivative_p_2 = sp.diff(substitute_m_into_coefficients_2_derivative_p_1, p)
substitute_m_into_coefficients_2_derivative_p_3 = sp.diff(substitute_m_into_coefficients_2_derivative_p_2, p)
substitute_m_into_coefficients_2_derivative_p_4 = sp.diff(substitute_m_into_coefficients_2_derivative_p_3, p)
#From the following equation, we know that the substitute_m_into_coefficients_2 is monotonically increasing with respect to p when p >= 1.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_p_1.subs(p, p_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_2 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_p_2.subs(p, p_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_3 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_p_3.subs(p, p_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_4 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_p_4.subs(p, p_lower_bound)),a)
print(f'\nThe first derivative function of substitute_m_into_coefficients_2 with respect to p: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_p_1),p)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_1: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_1}')
print(f'The second derivative function of substitute_m_into_coefficients_2 with respect to p: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_p_2),p)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_2: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_2}')
print(f'The third derivative function of substitute_m_into_coefficients_2 with respect to p: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_p_3),p)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_3: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_3}')
print(f'The fourth derivative function of substitute_m_into_coefficients_2 with respect to p: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_p_4),p)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_4: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_p_4}')

#--------------------------------------------------------------------------------------------------------------------
#Substitute p_lower_bound = 1 into substitute_m_into_diff_psi_squared.
substitute_m_p_into_diff_psi_squared = substitute_m_into_diff_psi_squared.subs(p,p_lower_bound)
#The formula of substitute_m_p_into_diff_psi_squared given in paper
substitute_m_p_into_coefficients_1 = (
    (324*a + 108)*t**4 + (648*a**2 + 4104*a + 1416)*t**3 +
    (432*a**3 + 4068*a**2 + 11284*a + 3292)*t**2 +
    (120*a**4 + 1264*a**3 + 6084*a**2 + 11136*a + 2744)*t +
    12*a**5 + 120*a**4 + 800*a**3 + 2640*a**2 + 3664*a + 768
)

substitute_m_p_into_coefficients_2 = (
    (81*a**2 + 2646*a + 981)*t**4 +
    (108*a**3 + 3618*a**2 + 14964*a + 4806)*t**3 +
    (54*a**4 + 1638*a**3 + 11845*a**2 + 28430*a + 7901)*t**2 +
    (12*a**5 + 278*a**4 + 2934*a**3 + 12864*a**2 + 22444*a + 5408)*t +
    a**6 + 12*a**5 + 220*a**4 + 1400*a**3 + 4564*a**2 + 6336*a + 1328
)

substitute_m_p_into_diff_psi_squared_in_paper = (substitute_m_p_into_coefficients_1*sp.sqrt(3*t + 3)+substitute_m_p_into_coefficients_2)
result_4 = sp.simplify(substitute_m_p_into_diff_psi_squared - substitute_m_p_into_diff_psi_squared_in_paper)
if result_4 == 0:
    print("\nIs the formula of substitute_m_p_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_p_into_diff_psi_squared in paper correct?: error")
# The above result will show that substitute_m_p_into_diff_psi_squared = substitute_m_p_into_diff_psi_squared_in_paper
# Both substitute_m_p_into_coefficients_1 and substitute_m_p_into_coefficients_2 can be regarded as functions of t.
# We need to show that both are monotonically increasing with respect to t.
# It is evident that substitute_m_p_into_coefficients_1 and substitute_m_p_into_coefficients_2 are quadratic functions of t.
# If the first, second and third derivatives of these two polynomials with respect to t are all greater than 0,
# then substitute_m_r_into_coefficients_1 and substitute_m_r_into_coefficients_2 will be monotonically increasing with respect to t.
substitute_m_p_into_coefficients_1_derivative_t_1 = sp.diff(substitute_m_p_into_coefficients_1, t)
substitute_m_p_into_coefficients_1_derivative_t_2 = sp.diff(substitute_m_p_into_coefficients_1_derivative_t_1, t)
substitute_m_p_into_coefficients_1_derivative_t_3 = sp.diff(substitute_m_p_into_coefficients_1_derivative_t_2, t)
#From the following equation, we know that the substitute_m_p_into_coefficients_1 is monotonically increasing with respect to t when t >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_1 = sp.collect(sp.simplify(substitute_m_p_into_coefficients_1_derivative_t_1.subs(t, t_lower_bound)),a)
substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_2 = sp.collect(sp.simplify(substitute_m_p_into_coefficients_1_derivative_t_2.subs(t, t_lower_bound)),a)
substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_3 = sp.collect(sp.simplify(substitute_m_p_into_coefficients_1_derivative_t_3.subs(t, t_lower_bound)),a)
print(f'The first derivative function of substitute_m_p_into_coefficients_1 with respect to t: {sp.collect(sp.expand(substitute_m_p_into_coefficients_1_derivative_t_1),t)}')
print(f'The value of substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_1: {substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_1}')
print(f'The second derivative function of substitute_m_p_into_coefficients_1 with respect to t: {sp.collect(sp.expand(substitute_m_p_into_coefficients_1_derivative_t_2),t)}')
print(f'The value of substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_2: {substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_2}')
print(f'The third derivative function of substitute_m_p_into_coefficients_1 with respect to t: {sp.collect(sp.expand(substitute_m_p_into_coefficients_1_derivative_t_3),t)}')
print(f'The value of substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_3: {substitute_lb_into_substitute_m_p_into_coefficients_1_derivative_t_3}')

substitute_m_p_into_coefficients_2_derivative_t_1 = sp.diff(substitute_m_p_into_coefficients_2, t)
substitute_m_p_into_coefficients_2_derivative_t_2 = sp.diff(substitute_m_p_into_coefficients_2_derivative_t_1, t)
substitute_m_p_into_coefficients_2_derivative_t_3 = sp.diff(substitute_m_p_into_coefficients_2_derivative_t_2, t)
#From the following equation, we know that the substitute_m_r_into_coefficients_2 is monotonically increasing with respect to t when t >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_p_into_coefficients_2_derivative_t_1 = sp.collect(sp.simplify(substitute_m_p_into_coefficients_2_derivative_t_1.subs(t, t_lower_bound)),a)
substitute_lb_into_substitute_m_p_into_coefficients_2_derivative_t_2 = sp.collect(sp.simplify(substitute_m_p_into_coefficients_2_derivative_t_2.subs(t, t_lower_bound)),a)
substitute_lb_into_substitute_m_p_into_coefficients_2_derivative_t_3 = sp.collect(sp.simplify(substitute_m_p_into_coefficients_2_derivative_t_3.subs(t, t_lower_bound)),a)
print(f'\nThe first derivative function of substitute_m_p_into_coefficients_2 with respect to t: {sp.collect(sp.expand(substitute_m_p_into_coefficients_2_derivative_t_1),t)}')
print(f'The value of substitute_lb_into_substitute_m_p_into_coefficients_2_derivative_t_1: {substitute_lb_into_substitute_m_p_into_coefficients_2_derivative_t_1}')
print(f'The second derivative function of substitute_m_p_into_coefficients_2 with respect to t: {sp.collect(sp.expand(substitute_m_p_into_coefficients_2_derivative_t_2),t)}')
print(f'The value of substitute_lb_into_substitute_m_p_into_coefficients_2_derivative_t_2: {substitute_lb_into_substitute_m_p_into_coefficients_2_derivative_t_2}')
print(f'The third derivative function of substitute_m_p_into_coefficients_2 with respect to t: {sp.collect(sp.expand(substitute_m_p_into_coefficients_2_derivative_t_3),t)}')
print(f'The value of substitute_lb_into_substitute_m_p_into_coefficients_2_derivative_t_3: {substitute_lb_into_substitute_m_p_into_coefficients_2_derivative_t_3}')

#--------------------------------------------------------------------------------------------------------------
#Substitute t_lower_bound = 3 into substitute_m_p_into_diff_psi_squared.
substitute_m_p_t_into_diff_psi_squared = substitute_m_p_into_diff_psi_squared.subs(t,t_lower_bound)
#The formula of substitute_m_p_t_into_diff_psi_squared given in paper
substitute_m_p_t_into_diff_psi_squared_in_paper = (
    a**6 + (24*sp.sqrt(3) + 48)*a**5 + (960*sp.sqrt(3) + 1540)*a**4 +
    (16960*sp.sqrt(3) + 27860)*a**3 + (150000*sp.sqrt(3) + 254008)*a**2 +
    (551360*sp.sqrt(3) + 947892)*a + 171216*sp.sqrt(3) + 297884
)
result_5 = sp.simplify(substitute_m_p_t_into_diff_psi_squared - substitute_m_p_t_into_diff_psi_squared_in_paper)
if result_5 == 0:
    print("\nIs the formula of substitute_m_r_t_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_r_t_into_diff_psi_squared in paper correct?: error")