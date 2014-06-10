from sympy import *
#k11=k21 k22=k12  and angle a01=a02 

c,a0,k1,k2,kappa,Er,Ez,Et,K=symbols("c,a0,k1,k2,kappa,E_r,E_z,E_theta,K" )
lambda_t,lambda_z,lambda_r,I_4,I_6,I_1,P,Ig=symbols("lambda_theta,lambda_z,lambda_r,I_4,I_6,I_1,P,I_g")
dSz,ints=symbols('\Delta\sigma_z,\int_ri^re')
init_printing(pretty_print=True)
#strech ratio
lt=sqrt(1+2*Et)
lz=sqrt(1+2*Ez)
lr=1/((lz)*(lt))
lr_Er=sqrt(1+2*Er)

#invariants
I1=lr**2+lt**2+lz**2
I1_Er=lr_Er**2+lt**2+lz**2      #First invariant with Er
I4=lz**2*(sin(a0))**2+lt**2*(cos(a0))**2

#energy function 

W=(c/2)*(I1-3)+(k1/(2*k2))*(exp(k2*(kappa*I1+(1-3*kappa)*I4-1)**2)-1)+(k1/(2*k2))*(exp(k2*(kappa*I1+(1-3*kappa)*I4-1)**2)-1)
#strain energy function including Er
W_Er=(c/2)*(I1_Er-3)+(k1/(2*k2))*(exp(k2*(kappa*I1_Er+(1-3*kappa)*I4-1)**2)-1)+(k1/(2*k2))*(exp(k2*(kappa*I1_Er+(1-3*kappa)*I4-1)**2)-1)

#kirchhoff stress 

S_tt=W.diff(Et)
S_zz=W.diff(Ez)
S_rr=W_Er.diff(Er)

#second derivative of strain energy 

S_tt_zz=S_tt.diff(Ez)
S_zz_zz=S_zz.diff(Ez)
S_rr_zz=W_Er.diff(Ez)

#simplify second drivatives 

S_rr_zz_stretch=((S_rr_zz.subs(2*Er+1,lambda_r**2)).subs(2*Et+1,lambda_t**2)).subs(2*Ez+1,lambda_z**2)
S_zz_zz_stretch=((S_zz_zz.subs(2*Ez+1,lambda_z**2)).subs(2*Et+1,lambda_t**2))
S_tt_zz_stretch=((S_tt_zz.subs(2*Ez+1,lambda_z**2)).subs(2*Et+1,lambda_t**2))

S_rr_zz_stretch_2=(S_rr_zz_stretch.subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_4)).subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_6).subs(lambda_r**2+lambda_z**+2+lambda_t**2,I_1)
S_zz_zz_stretch_2=(S_zz_zz_stretch.subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_4)).subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_6).subs((1/((lambda_z)*(lambda_t)))**2+lambda_z**+2+lambda_t**2,I_1)
S_tt_zz_stretch_2=(S_tt_zz_stretch.subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_4)).subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_6).subs((1/((lambda_z)*(lambda_t)))**2+lambda_z**+2+lambda_t**2,I_1)


#couchy Stress 

sigma_rr=lr**2*S_rr+K
sigma_zz=lz**2*S_zz+K
sigma_tt=lt**2*S_tt+K

#in terms of stretch ratios
sigma_rr_stretch=((sigma_rr.subs(2*Er+1,lambda_r**2)).subs(2*Et+1,lambda_t**2)).subs(2*Ez+1,lambda_z**2)
sigma_zz_stretch=((sigma_zz.subs(2*Ez+1,lambda_z**2)).subs(2*Et+1,lambda_t**2))
sigma_tt_stretch=((sigma_tt.subs(2*Ez+1,lambda_z**2)).subs(2*Et+1,lambda_t**2))



#replace I4 and I1 


sigma_rr_stretch_2=(sigma_rr_stretch.subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_4)).subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_6).subs(lambda_r**2+lambda_z**+2+lambda_t**2,I_1)
sigma_zz_stretch_2=(sigma_zz_stretch.subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_4)).subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_6).subs((1/((lambda_z)*(lambda_t)))**2+lambda_z**+2+lambda_t**2,I_1)
sigma_tt_stretch_2=(sigma_tt_stretch.subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_4)).subs(lambda_z**2*(sin(a0))**2+lambda_t**2*(cos(a0))**2,I_6).subs((1/((lambda_z)*(lambda_t)))**2+lambda_z**+2+lambda_t**2,I_1)


#sigma_r 
sigma_r=(sigma_tt_stretch_2-sigma_rr_stretch_2)
simple_sigma_r=simplify(((sigma_r.subs(lambda_t*lambda_z,1/lambda_r)).subs((I_1*kappa + I_4*(-3*kappa + 1) - 1),Ig)).factor(exp(Ig**2*k2)))

#sigma_z
sigma_z=sigma_zz_stretch_2-sigma_rr_stretch_2
simple_sigma_z=simplify(((sigma_z.subs(lambda_t*lambda_z,1/lambda_r)).subs((I_1*kappa + I_4*(-3*kappa + 1) - 1),Ig)).factor(exp(Ig**2*k2)))

#sigma_t 
sigma_t=sigma_tt_stretch_2-sigma_rr_stretch_2
simple_sigma_t=simplify(((sigma_t.subs(lambda_t*lambda_z,1/lambda_r)).subs((I_1*kappa + I_4*(-3*kappa + 1) - 1),Ig)).factor(exp(Ig**2*k2)))


#axial force 

N=2*sigma_zz_stretch_2-sigma_rr_stretch_2-sigma_tt_stretch_2
#simplify N replace (I_1*kappa + I_4*(-3*kappa + 1) - 1) with Ig
N_simple=simplify((simplify(N.subs((I_1*kappa + I_4*(-3*kappa + 1) - 1),Ig)).subs(lambda_t*lambda_z,1/lambda_r)).factor(exp(Ig**2*k2)))

#delta sigma_z 

del_sigma_z_1=simplify((2*(sigma_zz_stretch_2-K)*lambda_z**2*S_zz_zz_stretch_2-lambda_r**2*S_rr_zz_stretch_2).subs((I_1*kappa + I_4*(-3*kappa + 1) - 1),Ig))
simple_del_sigma_z_1=simplify(del_sigma_z_1.subs(lambda_t*lambda_z,1/lambda_r))

del_sigma_z_2=simplify((lambda_t**2*S_tt_zz_stretch_2-lambda_r**2*S_rr_zz_stretch_2).subs((I_1*kappa + I_4*(-3*kappa + 1) - 1),Ig))
simple_del_sigma_z_2=simplify(del_sigma_z_2.subs(lambda_t*lambda_z,1/lambda_r))


#del sigma_z replace del_sigma_z 

sigma_z_strain=(((((simple_sigma_z.subs(lambda_r**2,1+2*Er)).subs(lambda_z**2,1+2*Ez)).subs(lambda_t**2,1+2*Et)).subs(Ig,(I_1*kappa + I_4*(-3*kappa + 1) - 1))).subs(I_1,2*(Er+Et+Ez)+3)).subs(I_4,I4)
bbbb=-2*c*(2*Er + 1) + c*(2*Ez + 1) + 4*k1*(kappa*(2*Er + 2*Et + 2*Ez + 3) + (-3*kappa + 1)*((2*Et + 1)*cos(a0)**2 + (2*Ez + 1)*sin(a0)**2) - 1)*(-2*kappa*(2*Er + 1) - 3*kappa*(2*Ez + 1)*sin(a0)**2 + kappa*(2*Ez + 1) + (2*Ez + 1)*sin(a0)**2)*exp(k2*Ig**2)

#del_sigma_z_strain=simplify(((sigma_z_strain.subs(Er,0)).subs(Et,0)).subs(Ez,dSz))
del_sigma_z_strain=simplify(((bbbb.subs(Er,0)).subs(Et,0)).subs(Ez,dSz))
del_sigma_z_strain=simplify(del_sigma_z_strain.factor(dSz))

sigma_r_strain=(((((simple_sigma_r.subs(lambda_r**2,1+2*Er)).subs(lambda_z**2,1+2*Ez)).subs(lambda_t**2,1+2*Et)).subs(Ig,(I_1*kappa + I_4*(-3*kappa + 1) - 1))).subs(I_1,2*(Er+Et+Ez)+3)).subs(I_4,I4)
cccc=-2*c*(2*Er + 1) + c*(2*Et + 1) + 4*k1*(kappa*(2*Er + 2*Et + 2*Ez + 3) + (-3*kappa + 1)*((2*Et + 1)*cos(a0)**2 + (2*Ez + 1)*sin(a0)**2) - 1)*(-2*kappa*(2*Er + 1) - 3*kappa*(2*Et + 1)*cos(a0)**2 + kappa*(2*Et + 1) + (2*Et+ 1)*cos(a0)**2)
del_sigma_r_strain=simplify(((cccc.subs(Er,0)).subs(Et,0)).subs(Ez,dSz))
del_sigma_r_strain=simplify(del_sigma_r_strain.factor(dSz))

delta_sigma_z=del_sigma_z_strain+ints(ints(del_sigma_r_strain))