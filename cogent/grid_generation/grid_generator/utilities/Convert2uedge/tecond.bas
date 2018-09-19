# Transport coefficients
        kye = 1.0               #chi_e for radial elec energy diffusion
        kyi = 1.0               #chi_i for radial ion energy diffusion
        difni(1) = 1.0          #D for radial hydrogen diffusion
        travis(1) = 1.0         #eta_a for radial ion momentum diffusion
        difutm = 1.0            #toroidal diffusivity for potential

#-special values for the test
minu(1) = 2.    # ion mass relative to mp (hydrogen)
lnlam = 11.     # Coulomb logarithm
cthe = 0.       # thermal force coeff. for || mom. eq. (0.71 default)
cvgp = 0.       # turn off grad p in Te and Ti eqns.
cfvisx = 0.     # turn off viscous heating for ions
cfvisy = 0.     # turn off viscous heating for ions
ckinfl = 0.     # turn off viscous boundary term for heat flux
        
#-turn off convective heat flux 
cfcvte = 0.
cfcvti = 0.
cfloye = 0.
cfloyi = 0.

#turn off Y-flux coef for conv. in up-eq.
cmfy=0.0
isflxlde = 1
        
#turn off thermal force
cthe = 0.
        
#turn off equilibration
feqp = 0.

#turn off vgradP terms
cvgp = 0.

#turn off radiation energy gain/loss terms
chradi=0
chradr=0
