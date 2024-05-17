**************************************************************************
// Example I: Theoretical vs. empirical distribution
**************************************************************************

capture program drop sim_exp1
program define sim_exp1, rclass
	drop _all 
	local n = 1000
	local py_1 = .3
	local py_0 = .1
	local x_p = .2
	
	qui set obs `n'
	qui gen x = rbinomial(1, `x_p')
	qui gen y = .
	qui replace y = rbinomial(1, `py_1') if x == 1
	qui replace y = rbinomial(1, `py_0') if x == 0

	logit y x ,or
	return scalar est_b = _b[x]
	return scalar est_se_b = _se[x]
end

// simulate the study
sim_exp1 

// simulate the study 1,000 times
simulate est_b = r(est_b)  est_se_b = r(est_se_b) , seed(23016) reps(1000) : sim_exp1 

// summary statistics
summarize  est_b
di exp(r(mean))

* Expected mean and std deviation
di (60/140)/(80/720)
di sqrt(1/60+1/140+1/80+1/720)

scalar t_b = ln( (60/140)/(80/720) )
scalar t_se_b = sqrt(1/60+1/140+1/80+1/720)
di t_b+invnormal(.005)*t_se_b
di t_b+invnormal(.995)*t_se_b
di exp(t_b+invnormal(.005)*t_se_b)
di exp(t_b+invnormal(.995)*t_se_b)
su 

// figure
twoway ///
   (function t_b = normalden(x, t_b, t_se_b), range(`=t_b-4*t_se_b' `=t_b+4*t_se_b') lwidth(thick) color(black%80)) ///
   (hist est_b, bin(40) color(black%20)) ///
	, ytitle("Sampling Distribution") ylab(, nogrid) /// 
	xtitle("Estimated odds ratio (log scale)") ///
	legend(label(1 "Theoretical") label(2 "Simulated") ring(0) pos(1) col(1) region(style(none))) ///
	xlabel(`=ln(3.86)' "3.86" .85 "2.33" 1.85 "6.36", nogrid) plotregion(style(none)) graphregion(color(white)) ///
	aspect(0.8) name(figure_1, replace)	

***************************************************************************
// Example II: Confounding
***************************************************************************

// Introduce a confounding variable
clear all
*set seed 100723
set obs 1000
local p_y = .1
local p_x_z0 = .3
local p_w = .5
local p_z = .4

// w is an outcome predictor but not a confounder 
gen w = rbinomial(1, `p_w')

// z is a confounding variable because common cause of exposure and outcome 
gen z = rbinomial(1, `p_z')
gen x = rbinomial(1, invlogit(logit(`p_x_z0')+log(3)*z))
gen y = rbinomial(1, invlogit(logit(`p_y')+log(1.5)*x+log(2)*w+log(6)*z))

// Model 0. Omit w and z 
logit y x , or

// Model 1. Include an outcome predictor w but not the confounder z 
logit y x w, or

// Model 2. Include only the confounder z but not the outcome predictor w
logit y x z , or

// Model 3. Include both the confounder z and the outcome predictor w
logit y x z w , or

// simulating the unadjusted and adjusted effect for z a many times
capture program drop sim_exp2
program define sim_exp2, rclass
	drop _all 
	local n = 1000 
	local p_y = .1
	local p_x_z0 = .3
	local p_w = .5
	local p_z = .4
	set obs `n'
	gen w = rbinomial(1, `p_w')
	gen z = rbinomial(1, `p_z')
	gen x = rbinomial(1, invlogit(logit(`p_x_z0')+log(3)*z))
	gen y = rbinomial(1, invlogit(logit(`p_y')+log(1.5)*x+log(2)*w+log(6)*z))

	// Model 0. Omit w and z 
	logit y x , or
	ret scalar est_b1_m0 = _b[x]
	ret scalar est_se_b1_m0 = _se[x]
	
	// Model 1. Include an outcome predictor w but not the confounder z 
	logit y x w, or
	ret scalar est_b1_m1 = _b[x]
	ret scalar est_se_b1_m1 = _se[x]
	
	// Model 2. Include only the confounder z but not the outcome predictor w
	logit y x z , or
	ret scalar est_b1_m2 = _b[x]
	ret scalar est_se_b1_m2 = _se[x]

	// Model 3. Include both the confounder z and the outcome predictor w
	logit y x z w, or
	ret scalar est_b1_m3 = _b[x]
	ret scalar est_se_b1_m3 = _se[x]
	test x 
	ret scalar p_value_model3 = r(p)
end 

simulate ///
	est_b1_m0 = r(est_b1_m0) est_se_b1_m0 = r(est_se_b1_m0) ///
	est_b1_m1 = r(est_b1_m1) est_se_b1_m1 = r(est_se_b1_m1) ///
	est_b1_m2 = r(est_b1_m2) est_se_b1_m2 = r(est_se_b1_m2) ///
	est_b1_m3 = r(est_b1_m3) est_se_b1_m3 = r(est_se_b1_m3) ///
	p_value_model3 = r(p_value_model3) ///
	,  reps(1000): sim_exp2 // seed(140723)

// Power
count if p_value_model3 <0.05 
di r(N)/c(N)*100 

summarize  est_b1_m0, d
local mean_m0 = r(mean)
summarize  est_b1_m3, d
local mean_m3 = r(mean)

// figure 
twoway /// 
	(kdensity est_b1_m0, color(black) lpattern(solid) lwidth(medthick)) ///
	(kdensity est_b1_m1, color(black) lpattern(shortdash) lwidth(medthick)) ///
	(kdensity est_b1_m2, color(black) lpattern(longdash) lwidth(medthick)) ///
	(kdensity est_b1_m3, color(black) lpattern("-.-.") lwidth(medthick)), ///
	xline(`mean_m0' `mean_m3', lpattern(solid) lc(black%40)) ///
	plotregion(style(none)) graphregion(color(white)) /// 
	xlab(`mean_m0' `mean_m3', nogrid format(%3.2f) labsize(small) notick labgap(small)) /// 
	ytitle("Sampling distribution", size(small)) /// 
	xtitle("Estimated adjusted log odds ratio", size(small)) ///
	ylab(#8, angle(horizontal) nogrid notick labgap(small) labsize(small))  ///
	legend(label(1 "Model 0: Unadjusted") label(2 "Model 1: Outcome predictor but no confounder") label(3 "Model 2: Confounder but not outcome predictor") ///
	label(4 "Model 3: Confounder and outcome predictor") ring(1) col(1) pos(6)) /// 
	name(figure_2, replace) scheme(s1mono)
	
***************************************************************************
// Example III: Missing Data
***************************************************************************
clear all
* set seed 1151
set obs 1000
local p_y = .1
local p_x = .3
local p_z = .4

gen z = rbinomial(1, `p_z')
gen x = rbinomial(1, invlogit(logit(`p_x')+log(3)*z))
gen y = rbinomial(1, invlogit(logit(`p_y')+log(1.5)*x+log(6)*z))

logit y x z 
est store full

// MCAR for the confounder 70% 
replace z = . if runiform()<.7

logit y x z 
est store miss

// impute missing variable
mi set wide
mi register imputed z 
mi impute logit z x y, add(70)

mi estimate, post: logit y x z 
est store imp 

// compare full, complete case, and imputed analysis
est table full miss imp, b se 

// simulation 
capture program drop sim_exp3
program define sim_exp3, rclass
drop _all 
	
	set obs 1000 
	local p_y = .1
	local p_x = .3
	local p_z = .4
	gen z = rbinomial(1, `p_z')
	gen x = rbinomial(1, invlogit(logit(`p_x')+log(3)*z))
	gen y = rbinomial(1, invlogit(logit(`p_y')+log(1.5)*x+log(6)*z))
	
	logit y x z 
	ret scalar x_full = _b[x]
	ret scalar x_se_full = _se[x]
	
	replace z = . if runiform()<.7  
	logit y x z 
	ret scalar x_miss = _b[x]
	ret scalar x_se_miss = _se[x]
	
	mi set wide
	mi register imputed z 
	mi impute logit z x y, add(70)
	mi estimate, post: logit y x z 
	ret scalar x_imp = _b[x]
	ret scalar x_se_imp = _se[x]
end 

simulate x_full = r(x_full) x_se_full = r(x_se_full) x_miss = r(x_miss) x_se_miss =r(x_se_miss) x_imp = r(x_imp) x_se_imp = r(x_se_imp), reps(1000): sim_exp3

foreach m in x_full x_miss x_imp{
	su `m'
	local `m' = r(mean)
}

tw ///
	(kdensity x_full, color(black) lpattern(solid) lwidth(medthick)) ///
	(kdensity x_miss, color(black) lpattern(dash) lwidth(medthick)) ///
	(kdensity x_imp, color(black) lpattern(shortdash) lwidth(medthick)), ///
	plotregion(style(none)) graphregion(color(white)) ///
	legend(label (1 "Full analysis") label (2 "Complete case analysis") label (3 "Analysis after imputation") ///
	ring(1) col(1) pos(6) size(small))  /// 
	xtitle("Estimated adjusted log odds ratio", size(small)) ytitle("Sampling distribution", size(small)) /// 
	xlab(`x_full', nogrid format(%3.2f) labsize(small) notick labgap(small)) ///
	xline(`x_full', lpattern(solid) lc(black%40)) /// 
	ylab(#8, labsize(small) angle(horizontal) nogrid notick labgap(small)) /// 
	scheme(s1mono) name(figure_3, replace)
