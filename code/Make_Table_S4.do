
* This file creates Table S4 in Hsiang Burke Miguel 2013

clear all

* replace the following line with the directory where the replication data was unzipped
global wd "\Documents\Dropbox\Marshall-Sol\drafts\Science_review_v4\replication"

insheet using $wd/data/standardized_effects.csv, comma clear

* main calculation uses authors' estimates of t-stats, not ours
gen tstat= abs(their_effect/their_se)
gen logtstat = log(tstat)
gen logdof = log(sqrt(their_dof))
gen temp = (independentvariable=="ENSO"| independentvariable=="PSDI" | independentvariable=="SPEI" | independentvariable=="temperature")

reg logtstat logdof
reg logtstat logdof if temp==1


