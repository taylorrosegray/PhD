from object_library import all_form_factors, FormFactor
from function_library import complexconjugate, re, im, csc, sec, acsc, asec, theta_function

### form factor from arXiv: 0906.0580 eqs A18 + A19
G2 = FormFactor(name = 'G2',
                 type = 'real',
                 value = '(((111*74**(-1./3.)/Me)**2*(-P(-1,1)*P(-1,1)))/(1+(111*74**(-1./3.)/Me)**2*(-P(-1,2)*P(-1,2))))**2*(1/(1+(-P(-1,3)*P(-1,3))/(0.164*184.**(-2./3.)))**2)*74**2\
                  + ((773*74**(-2/3)/Me)**2*(-P(-1,4)*P(-1,4))/(1+(773*74**(-2/3)/Me)**2*(-P(-1,5)*P(-1,5))))**2 * ((1+((-P(-1,6)*P(-1,6))/(4*0.938**2)*(2.79**2 - 1)))/(1+(-P(-1,6)*P(-1,6))/0.71)**4)**2 * 74') # elastic + inelastic G2
