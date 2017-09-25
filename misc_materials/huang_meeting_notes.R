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