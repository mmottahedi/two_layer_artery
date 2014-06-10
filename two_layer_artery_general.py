from sympy import *
#general case

c,a01,a02,k11,k12,k21,k22,kappa,Er,Ez,Et,K=symbols("c,a01,a02,k11,k12,k21,k22,kappa,Er,Ez,Et,K" )
lambda_t,lambda_z,lambda_r,I_4,I_6,I_1,P=symbols("lambda_theta,lambda_z,lambda_r,I_4,I_6,I_1,P")
init_printing(pretty_print=True)
#strech ratio
lt=sqrt(1+2*Et)
lz=sqrt(1+2*Ez)
lr=1/((lz)*(lt))
lr_Er=sqrt(1+2*Er)

#invariants
I1=lr**2+lt**2+lz**2
I1_Er=lr_Er**2+lt**2+lz**2      #First invariant with Er
I4=lz**2*(sin(a01))**2+lt**2*(cos(a01))**2
I6=lz**2*(sin(a02))**2+lt**2*(cos(a02))**2

#energy function 

W=(c/2)*(I1-3)+(k11/(2*k12))*(exp(k21*(kappa*I1+(1-3*kappa)*I4-1)**2)-1)+(k21/(2*k22))*(exp(k22*(kappa*I1+(1-3*kappa)*I6-1)**2)-1)
#strain energy function including Er
W_Er=(c/2)*(I1_Er-3)+(k11/(2*k12))*(exp(k21*(kappa*I1_Er+(1-3*kappa)*I4-1)**2)-1)+(k21/(2*k22))*(exp(k22*(kappa*I1_Er+(1-3*kappa)*I6-1)**2)-1)

#kirchhoff stress 

S_tt=W.diff(Et)
S_zz=W.diff(Ez)
S_rr=W_Er.diff(Er)

#couchy Stress 

sigma_rr=lr**2*S_rr+K
sigma_zz=lz**2*S_zz+K
sigma_tt=lt**2*S_tt+K

#in terms of stretch ratios
sigma_rr_stretch=((sigma_rr.subs(2*Er+1,lambda_r**2)).subs(2*Et+1,lambda_t**2)).subs(2*Ez+1,lambda_z**2)
sigma_zz_stretch=((sigma_zz.subs(2*Ez+1,lambda_z**2)).subs(2*Et+1,lambda_t**2))
sigma_tt_stretch=((sigma_tt.subs(2*Ez+1,lambda_z**2)).subs(2*Et+1,lambda_t**2))

#replace I4 and I6 

sigma_rr_stretch_2=(sigma_rr_stretch.subs(lambda_z**2*(sin(a01))**2+lambda_t**2*(cos(a01))**2,I_4)).subs(lambda_z**2*(sin(a02))**2+lambda_t**2*(cos(a02))**2,I_6).subs(lambda_r**2+lambda_z**+2+lambda_t**2,I_1)
sigma_zz_stretch_2=(sigma_zz_stretch.subs(lambda_z**2*(sin(a01))**2+lambda_t**2*(cos(a01))**2,I_4)).subs(lambda_z**2*(sin(a02))**2+lambda_t**2*(cos(a02))**2,I_6).subs((1/((lambda_z)*(lambda_t)))**2+lambda_z**+2+lambda_t**2,I_1)
sigma_tt_stretch_2=(sigma_tt_stretch.subs(lambda_z**2*(sin(a01))**2+lambda_t**2*(cos(a01))**2,I_4)).subs(lambda_z**2*(sin(a02))**2+lambda_t**2*(cos(a02))**2,I_6).subs((1/((lambda_z)*(lambda_t)))**2+lambda_z**+2+lambda_t**2,I_1)


#sigma_tt-sigma_rr
sigma_t_min_sigma_r=(sigma_tt_stretch_2-sigma_rr_stretch_2)
#sigma_t_min_sigma_r=simplify(sigma_tt_stretch_2-sigma_rr_stretch_2)

#axial force 

N=2*sigma_zz_stretch_2-sigma_rr_stretch_2-sigma_tt_stretch_2