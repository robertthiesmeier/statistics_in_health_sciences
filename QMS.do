// Stata file for QMS Monday 20 May 2024

***************************************************************
// coin flip 
***************************************************************
local c = rbinomial(1,0.5)
di `c'
di in red cond(`c'<1, "Heads", "Tails")

// many coin flips
local reps = 10
local t = 0 
local i = 1
while `i' < `reps'{
	local c = rbinomial(1,0.5)
	if `c'<1 local t = `t'+1
	local i = `i' + 1
}

di in red "Probability of Tails = " `t'/`reps'

cap prog drop coin_flip
prog define coin_flip, rclass 
	drop _all 
	set obs 100
	gen c = rbinomial(1,0.5)
	su c 
	ret scalar flip = r(mean)
end 
coin_flip
simulate flip = r(flip), reps(1000): coin_flip
hist flip , discrete xtitle("Coin flip") color(orange%50) xline(0.5) xlabel(0.4 0.5 0.6)
graph export coin_flip.png, replace width(4000)


***************************************************************
// 2 by 2 table
***************************************************************
cap program drop or 
program define or , rclass 
	drop _all
	set obs 1000
	gen x = rbinomial(1, 0.2)
	gen y = rbinomial(1, 0.1 + (0.3-0.1)*x)
	cs y x , or
	ret scalar or = r(or)
end 

// single table
or
            
simulate or = r(or), reps(1) seed(19052024): or

// many replications
simulate or = r(or), reps(1000): or

hist or , xtitle("Odds Ratio") color(blue%50) xline(3.9) xlabel(3.9)

***************************************************************
// simple confounding mechanism
***************************************************************
cap prog drop confounding
program define confounding, rclass 
	drop _all 
	
	set obs 1000
	
	gen c = rbinomial(1, 0.3)
	gen x = rbinomial(1, invlogit(logit(0.4)+ln(3)*c))
	gen y = rbinomial(1, invlogit(logit(0.05)+ln(1.8)*x+ln(5)*c))
	
	// undadjusted
	logit y x, or
	ret scalar b1_unadjusted = _b[x]
	
	// adjusted effect 
	logit y x c, or 
	ret scalar b1_adjusted = _b[x]
	test x 
	ret scalar p = r(p)
	
end 

// single analysis
*set seed 1
*confounding

// simulate many studies
simulate b1_unadjusted = r(b1_unadjusted) b1_adjusted = r(b1_adjusted) p = r(p), reps(1000): confounding 

count if p < 0.05 
di "Power to detect an adjusted effect = " r(N)/c(N) * 100 "%"

gen or0 = exp(b1_unadjusted)
gen or1 = exp(b1_adjusted)

qui su or0
local mean_or0 = r(mean)
qui su or1
local mean_or1 = r(mean)

twoway /// 
	(kdensity or0, recast(area) color(orange%40)) /// 
	(kdensity or1, recast(area) color(blue%40)), /// 
	xtitle("Odds ratio (log)") xlabel(`mean_or0' `mean_or1', format(%3.2f)) /// 
	xline(`mean_or0' `mean_or1') xscale(log) /// 
	legend(label(1 "Unadjusted") label(2 "Adjusted") pos(6) row(1)) name(fig1, replace)
