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
## a. maybe try using a binomial or poisson model on the data - NOTEDs
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
##     he can't upload to subdirectories
## b. Ray will move virtual pod to sachsURAP
## c. Make the 3.2.4 CI overlay plot
## d. Plot CI overlay without calibrated coefficients
## e. If need help -> Yimin -> Peter -> Andy in that order
## f. Downsize plot legends, assume audience will view only in IDE
## 
