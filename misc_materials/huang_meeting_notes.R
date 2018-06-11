## SACHS, EBERT & HUANG @ STRADA 09.06.17 
# 1. write script to calculate I(d) using incremental effect additivity for a 
#    low LET IDER and N >= 1 HZE IDERs
#    a. stay below 10 gray, but look beyond 1 gray for mathermatical reasons
#    b. look at 6-7 gray for low LET IDERs
#    c. 0-1 gray for high LET ions (HZE)

## In my minutes of Wed. 9/6/17 meeting I forgot to add an additional assignment
## that Edward and I agreed on. At the moment merge2 has a function to calculate
## baselines for mixtures of any number N of HZE and one to calculate baselines 
## for a mixture of one HZE with one low LET ion. The latter should be extended 
## to mixtures with N HZE ions and one low LET ions. Maybe we only need one R 
## function to calculate simple and incremental effect additivity baselines for 
## N>=1 HZE and either 0 or 1 low LET ion.

## SACHS & HUANG @ STRADA 09.14.17
## Read methods, results of Dae paper

## SACHS & HUANG @ STRADA 09.21.17
## debug yimin's code 
## a. runtime of yimins code
## b. original merge2 buggy?
## make new file with merge2 called lowhighEdoseExploration - DONE
## a. only give ray permission if possible 
## look into physics for politicans if ebeling doesn't add - DONE
## provide instructions for making rStudio project with gitHub on the issues page
## Sys.time() on issues page - DONE
## apply vs for loop on issues page - DONE

## SACHS & HUANG @ EVANS 02.23.18
## 7 High LET and HZE plots (4 top 3 bottom) and one low HZE plot bottom right (single ion beams)
## Lori Chappelle (Cucinotta's former statistician) wants to collaborate (She'll meet with Ray on Wednesday)
## I will hold off on edits until after Yimin's latest updates are live on GitHub
## Change HZE and LET acronyms in cript to upper-case - DONE

## SACHS & HUANG @ EVANS 03.02.18
## Will submit to REB instead of LSSR
## Next batch of mice will be irradiated with "representative beams"
## Try to get close to the millimeter figure guidelines
## Won't need new LETs for the minor paper
## Take the path of least resistance when making figures (either use existing code or overhaul first)
## We're not worrying about TE only models
## We're not worrying about older models
## Need to rewrite calculate I(d) for clarity

## 03.15.18 SACHS & HUANG 918 EVANS
## Last night Bragg peaks calculations 
## i. May trump biophysical reactions
## ii. As ions slow down nuclear collisions become more important
## iii. LET increases exponentially as ion slows
## iv. Current radiation studies may be misguided
## v. Correlation if points near the end of their path show abnormally large tumorgenesis
## How did the meetings go?
## Discuss new materials
## Will need to fix calculate_i(d) later to take arbitrary number of low-LET ions and dose
## Need to also fix calculate_SEA to take low-LET IDERs - DONE
## This current paper excuses us from having to do MIXDERs with over two low-LET IDERS
## Communicate through master file (i.e. comments on .plottingYimin) - NOTED
## Stick to .R and base plot - NOTED
## Commit more often - NOTED
## Make some IDER lines dashed to show IDERs covered underneath - DONE
## Fix calculate_SEA to be correct - DONE
## Plots will change over time - NOTED
## Hold off on editing Bragg script for now - NOTED
## No more emailing .eps - stick with having a coding block. - NOTED
## 4-panel - NOTED


## SACHS & HUANG 918 EVANS 032318
## size all legends to native plot viewer - DONE
## make sure to keep repo permissions clear at the organizaton level - NOTED
## write Ray an email reminder to ask Yimin to refactor his Monte Carlo code from Feb. - DONE
## ask Yimin to refactor his Monte Carlo code  - NOTED
## error bars for the 8panel, 2 * standard deviation or 95% confidence intervals - NOTED
## a. maybe try using a binomial or poisson model on the data - NOTED
## get n x n normalized correlation matrix from vcov() matrix  - DONE
## potential big project: refactor Andy's CA code with Peter 
## error bars for plots - NOTED
## launch github organization - IN PROGRESS

## SACHS & HUANG @ THE STRADA 040118
## Sachs checked NWEIGHT, values are all wrong, Sachs has to change - NOTED
## If possible maybe make monteCarlo.R more efficient 
## Produce ribbon plots (with ggplot2) of Monte Carlo in plots.R - DONE
## May have to eventually write special script or add to current HG script
## Look into finding out how to make data points translucent - DONE
## Sachs will figure out if photon model is correct - NOTED
## Long run should use Bayesian or Maximum-likelihood instead of nls()
## NOTES ON LET AND EFFECTIVENESS: - NOTED
## i. Must take into account stochastic energy deposition of ions
## ii. For example low-LET ions can be described with a simple Poisson model
##     while li-LET ions are more similar to a complex Poisson model which 
##     results in more effectiveness

## Foreseeable long run tasks:
## a. Have more than one lowLET IDER (possibly characterized by LET, maybe other parameter)
## b. Need TE-only models and figures
## c. Compare NTE-TE and TE-only models, use some kind of metric (IC, machine 
##    learning, goodness of fit, etc.)
## d. Have to allow for advanced Incremental-Effect Additivity using autonomous 
##    differential equations and one-ion models characterized by their slopes as
##    functions of effect.

## 040618 SACHS & HUANG 918 EVANS
## Ignore Beta Decay and Control (dose = 0) data
## From now on data will be organized into ion_data, beta_decay, and control
## Ray will send his latest local copy synergyTheory.R 
## i. Incorporate changes to GitHub master branch
## Script approaching final versions
## Main bottleneck now is paper draft
## Ray will send Web Supplement, writing is optional

## 041318 SACHS & HUANG 918 EVANS / YALIS
## a. Move MS06 to paper, check Ray's permissions b/c 
##     he can't upload to subdirectories - CHANGED PERMISSIONS
## b. Ray will move virtual pod to sachsURAP - NOTED
## c. Make the 3.2.4 CI overlay plot - DONE
## d. Plot CI overlay without calibrated coefficients - DONE
## e. If need help -> Yimin -> Peter -> Andy in that order - NOTED
## f. Downsize plot legends, assume audience will view only in IDE - DONE

## 042018 SACHS & HUANG 918 EVANS / YALIS
## Paper is getting shorter; progress is slow; - NOTED
## Consider the parameter order dose, LET, z, ratios (z = lowLET) - DONE
## Change calculate_complex_id to calculate_id - DONE
## Test calculate_id on empty MIXDER - DONE, error in uniroot
## Change name of hgData.R to dataandComments or equivalent - DONE
## consider maybe longer commenting blocks in synergyTheory - NOTED
## Consider changing "corr" in monteCarlo.R to respect_corr or something equivalent try "vcov" - DONE
## look at URAP.CA monte carlo in "Graphs" 

## 042718 SACHS & HUANG 918 EVANS 
## Add summary call to line 54 of plots.R. This would be on the low LET model coefficient - DONE
## Change calib_ in DER names to calibrated_ - DONE
## See if I can find any examples of (dose, ratios, LET, ...) in function headers - DONE
## All Fe56 600 MeV/u LET values at 195 need to go back to 193 - DONE
## Use yellow for vcov confidence interval (inside one) - DONE
## Use aquamarine (aquamarine2), lightblue, for outside non-vcov confidence interval - DONE
## Look into getting Illustrator - NO LONGER RELEVANT
## Look into getting mathtype onto linux - NO LONGER RELEVANT

## 050418 SACHS & HUANG STRADA CAFE
## Ray will try using LaTeX and Illustrator - NOTED
## Edward will complete past open assignments - DONE
## Iff Claire is interested in EvoLab then Edward will connect her with current
## lab members. - NOTED


## 0501618 SACHS & HUANG GMAIL
# 1.
# Remove all figures not in the newest version of the minor paper (ms09.docx)
# or its online supplement from plot.R - DONE BY RAY

# 2.
# Consider defining dose vectors such as forty_cGy etc. with each figure 
# instead of in the global environment. With only about 6 figures using this 
# notation, there will be only minor redundancies. The advantage would be that
# more essential information would be available right at the figure and could be
# tweaked by changing the dose vector in a way unique to to optimizing that 
# figure, e.g. by adding one extra very low dose to the dose vector as suggested
# below. - DONE BY RAY

# 3.
# In plots.R renumber Figs. 1-11 as in the attached paper (some are not in 
# plots.R but rather in minor auxiliary .R files; just skip those and go 
# directly to the next Fig. that is in plots.R)  - DONE BY RAY

# 4.
# Order figures consequtively within plots.R. Perhaps put online supplement 
# figs at end or in a separate file? - DONE BY RAY

# 5.
# Many of the plots have the slope at the origin much too small. If a dose 
# axis runs from 0 to 50 or 100 Gy the region near 0 should look as if the
# slope at zero is infinite (despite the fact that it is actually finite but 
# very large) and should appear to have at discontinuous first derivative (which
# is actually a region of very high concavity) at a kink. These misimpressions 
# will occur if the first non-zero dose point is, say, 1 cGy. There are various 
# reasons why one might want extra dose points near the origin. If there is 
# no other reason for extra dose points adding a single dose point at, say 
# 10^-5 Gy (i.e.10^-3 cGy, a thousand times smaller) is enough to cure the 
# misimpression. Just make the IDER effect or mixture DER at the extra dose 
# point a little bit smaller than the IDER or mixture effect calculated near 
# 1 cGy (after the usual subtracting out of background). Then the slope will 
# look as if it might be infinite and the kink will look as if it could be a 
# slope discontinuity as desired. No further calculations are then required. 
# This trick can and should be used iff: (a) the maximum dose on the graph is 
# larger than a few cGy; and (b), there is no other reason (such as feeding 
# uniroot) to have lots of extra dose points near the origin. - NOTED

# 6.
# At the moment the command "lowLET='True'" need lots of attention.
# It cannot be used if one has both protons and alpha particles because there 
# is only one DER for both protons and alpha particles and the alpha particles 
# would be incorrectly assigned an HZE IDER with LET 1.6. If we are calculating
# I(d) for a mixture that has only HZE components then the opposite problem 
# occurs: the first HZE component will (correct me if I am wrong) be incorrectly
# assigned the low LET IDER and I(d) will be somewhat too small near dose=0. I 
# think the two ribbon plots (now Figs. 10 and 11) are both wrong, not only in 
# the I(d) curve but also in the ribbon widths, for this reason. - NEED TO 
# DISCUSS WITH RAY

# 7. 
# In many places "calibrated" is written where "correlated" is meant.
# (plots.R, check all monte_carlo_runs) - DONE

# 8.
# In at least one place, the "corrTop" and "corrBottom" curves have been 
# switched in my sandbox, so that the former is lower than the latter instead 
# of being above it. Please check whether the master version on gitHub of 
# plot.R also has this problem (error is in line 222, bug probably in 
# monte_carlo.R) - FIXED

# 9.
# The maximum doses in Figs. Now numbered 10 and 11 are too large. The Figs 
# should go up to 40 or 50 cGy, the same as their non-ribbon counterparts in 
# section 3.2 of ms09 (Fig. 8 and Fig. 9 IIRC). - DONE BY RAY 

# 10.
# I will use plot() without opacity commands because the opacity commands are 
# giving me trouble. Fig. 11 in MS09.docx shows the kind of  ribbon figure we 
# need IMHO but never used layers or opacity. Instead the order of the various
# curves was chosen such that the later give the desired effects when they 
# overwrite the earlier curves. You can make the master with ggplot2() and/or 
# layers and/or using opacity commands if you want but the results for ribbon 
# plots should look similar to Fig. 11 in MS09.docx. - NOTED 

# 11. (Additional Bug)
# I just saw that when we went over to the .csv way of handling data we 
# reordered the data using Z as a label (first protons with Z=1, then Helium 
# .with Z=2, then neon with Z=10, then Silicon with Z=14,  and so on up. But we
# did not correct synergy.R and plot.R for this reordering (I don't know if 
# it matters to Monte_Carlo.R?). As far as I can tell at the moment, 
# the only change is for the low LET model: 
# on our graphs all of the proton points are mislabelled as helium and four of
# the helium points are mislabeled as protons. Otherwise we were often smart 
# enough to use names instead of numbers so hopefully the mislabelling is the
# only mistake, all the regressions and calculations of I(d) are probably OK 
# (though the latter is in trouble when low_LET=TRUE is misused, as I wrote 
# you earlier.) - CORRECTED BY RAY


## 053118 SACHS & HUANG @ THE STRADA 
# Make 8 panel plot of HZE ions with points and error bars,
# follow procedure from Ray's Figure 5, base R or ggplot no preference


# Add low-LET boolean parameter in simulate_monte_carlo
# (only if below is not implemented) - NEGATED BY BELOW

# MODIIFY calculate_id
# For a MIXDER of N low-LET ions and 0 HZE ions, add up the LETs of all 
# the low-LET ions and treat as IDER (Z < 3 is low-LET)
# For a MIXDER of M low-LET ions and N HZE ions, add up the LETs of all 
# the low-LET ions and treat as MIXDER of N HZE + one low-LET DER (Z < 3 is low-LET)
# ( IF DONE CORRECTLY THERE IS NO LONGER A NEED FOR lowLET boolean toggle parameters)
# - DONE

# If par() is ever called please reset the graphics settings - DONE

# If extra time rewrite Mark's report and code on Lam's synergy theory - LATER

# Ray will send me link to audrey's repo about the giant letters - LATER


