# Correcting Alpen data for LET change between beam injection and mouse skin. 
#Under construction but probably bug-free. Here tracks in water of almost fully ionized 2He4 atoms (alpha particles) are analyzed as an example.
# Many physics approximations are used in the calculations but the main source of uncertainty is biological variability, not these approximations.
# The main quantitites of interest are kinetic energy in MeV/u=ke, LET in keV/u, range=average straight line distance travelled before an alpha particle stops at its Bragg peak. Measured in cm
rm(list=ls()) #The following data is taken from NASA's GERM software, a program that can calculate properties of radiation tracks. We do not need the full power of the GEANT4 suite for our estimates.
he4_ke_x=c(60,75,90,100,125,150,175,200,250,300,350) # This is kinetic energy ke; but smooth.spline wants it to be "x" and not any other character string. So he4_ke_x will soon be replaced by x
he4_range=c(3.07,4.58,6.34,7.65,11.37,15.65,20.48,25.78,37.68,51.09,65.88) # The corresponding ranges
he4_LET=c(4.35,3.66,3.18,2.94,2.49,2.19,1.97,1.81,1.57,1.42,1.30) #The corresponding LETs
# JE Turner Atoms, Radiation, and Radiation Protection. 3d edition. Eq. 5-42 says that if 2 ions have equal MeV/u then to approximation adequate for our purposes (range_1/range_2)=(M_1/M_2)*(Z_2/Z_1)^2
# e.g.(range_1/range_Fe56)=(M_1/56)*(26/Z_1)^2 if ion2 is 26Fe56.
# Evntually we will get an estimate of how much matter was between Alpen's upstream ion injections and mouse skin from comparing Alpen's and Chang's data on 26Fe56 injected with about 600 MeV/u.
# Then we can use that estimate to correct Alpen's LETs at the mouse skin for Alpen's proton, He4,Ne20, and Fe56 at 350 MeV/u data by using Turner's Eq.5-42.
x=he4_ke_x 
range_of_ke=smooth.spline(x,he4_range,df=5) # list of 15 #df must be less than length(x)=length(range)# list of 15
range_fctn_of_ke=function(ke)unname(predict(range_of_ke, x=ke)[[2]]) #extract function from list
LET_of_ke=smooth.spline(x,he4_LET,df=5) # similar to 2 lines above
LET_fctn_of_ke=function(ke)unname(predict(LET_of_ke, x=ke)[[2]])
x=he4_range 
ke=he4_ke_x
ke_of_range=smooth.spline(x,ke,df=10) # compositional inverse list of 15
ke_fctn_of_range=function(he4_range)unname(predict(ke_of_range, x=he4_range)[[2]]) # compositional inverse function
range=40:1
ke=ke_fctn_of_range(range)
LET=LET_fctn_of_ke(ke)
plot(40-range,LET,type='l') # This curve is merely the famous Bragg peak curve, and is merely an interpolation for the data we started with. But the smooth interpolation and the formula in Turner can give us the desired estimates of how much LET Alpen's ions gained before they hit the mice.