***********************************************
* Author: Tyson Franke
* Creation: Nov 19, 2025
* Summary: This is a do-file for the final project of econ 460
* Data Links: 
* https://fred.stlouisfed.org/series/LNS13023622
* https://fred.stlouisfed.org/series/USRECM
* https://fred.stlouisfed.org/series/AAA
* https://fred.stlouisfed.org/series/BAA
* https://fred.stlouisfed.org/series/GS1
* https://fred.stlouisfed.org/series/GS10
* https://fred.stlouisfed.org/series/TB3MS

* Input files: LNS13023622.csv, USRECM.csv, AAA.csv, BAA.csv, GS1.csv, GS10.csv, TB3MS.csv

* Temp files: 

* Output files: forecast.png, forecast_constant, forecast_analysis.png
* job_losers_sa.dta, USRECM.dta, AAA.dta, BAA.dta, GS1.dta, GS10.dta, TB3MS.dta
***********************************************

clear all
set more off
capture log close
global proj_dir "your-project-directory"
log using "${proj_dir}\forecast_analysis.log", replace text

* Setup
import delimited "${proj_dir}\LNS13023622.csv", clear
rename lns13023622 job_losers_sa
save "job_losers_sa.dta", replace
import delimited "${proj_dir}/USRECM.csv", clear
save "USRECM.dta", replace
import delimited "${proj_dir}/AAA.csv", clear
save "AAA.dta", replace
import delimited "${proj_dir}/BAA.csv", clear
save "BAA.dta", replace
import delimited "${proj_dir}/GS1.csv", clear
save "GS1.dta", replace
import delimited "${proj_dir}/GS10.csv", clear
save "GS10.dta", replace
import delimited "${proj_dir}/TB3MS.csv", clear
save "TB3MS.dta", replace

use job_losers_sa.dta, clear
local files : dir . files "*.dta"
foreach file  in `files' {
    merge 1:1 observation_date using `file', nogen
}
drop if missing(job_losers_sa)

gen t = mofd(date(observation_date, "YMD"))
format t %tm
tsset t
drop if t >= ym(2025, 9)

gen covid_shock = (t >= tm(2020m2) & t <= tm(2020m4))
gen month = month(dofm(t))

gen spread1 = gs1 - tb3ms
gen spread2 = gs10 -tb3ms
gen corporate = baa - aaa
gen dt3 = tb3ms - L.tb3ms
gen dt12 = gs1 - L.gs1

*--------------------
describe
list in 1/10
tsline job_losers_sa usrecm
*--------------------

/*
	Modeling
_______________*/
* Denote: SS = Statistically Significant, 
* SI = Statistically Insignificant

* Stationarity
dfuller job_losers_sa, regress // SS, _cons SI 
ac job_losers_sa // Slow Decay
pac job_losers_sa // Cutoff after lag 1
// Conflicting results could be due to SA

** No Unit Root Assumption
* Mean Shift
reg job_losers_sa covid_shock // SS
predict mean_shift, xb
predict demeaned_sa, residual
* Trend
reg demeaned_sa t // SS but t coef 0.0133067  
predict trend, xb
tsline trend job_losers demeaned_sa
predict demeaned_detrended_sa, residual
* Remaining Seasonality
reg demeaned_detrended i.month // SI
predict seasonal, xb
tsline seasonal // Marginal differences
* Stationarity || Demeaned & Detrended
dfuller demeaned_detrended_sa, regress // SS, _cons SI
ac demeaned_detrended_sa // Slow Decay
pac demeaned_detrended // Cutoff after 3
// Suggests AR(3)
arima demeaned_detrended_sa, arima(3,0,0)
predict resid_ar3, residual
ac resid_ar3 // White Noise
pac resid_ar3 // White Noise

* AR(3) + Trend + Recession/Covid Dummies
arima job_losers usrecm covid_shock t, arima(3,0,0) // t, usrecm are SI
* AR(3) + Trend + Covid Dummy
arima job_losers covid_shock t, arima(3,0,0) // t is SI
* AR(P) + Covid Dummy
forvalues i = 1/3 {
	arima job_losers covid_shock, arima(`i',0,0)
	estimates store AR`i'_model
}
arima job_losers covid_shock, arima(3,0,0)
predict resid_ar3_cov, residual
ac resid_ar3_cov // White Noise
pac resid_ar3_cov // White Noise

** Existing Unit Root Assumption
* Stationarity || Differencing
dfuller D.job_losers, regress // SS, _cons SI 
ac D.job_losers // White Noise
pac D.job_losers // White Noise
// Suggests Random Walk, I(1)


* Differencing Models
arima D.job_losers_sa
estimates store diff_model
arima D.job_losers_sa, noconstant
estimates store diff_model_no_cons


** Leading Indicators
dfuller spread1 // SS
dfuller spread2 // SS
dfuller corporate // SS
dfuller dt3 // SS
dfuller dt12 // SS

local indicators dt3 dt12 spread1 spread2 corporate

forvalues l = 1/20 {
    foreach var of local indicators {
		// Granger Causality
        quietly regress D.job_losers_sa L(1/`l').D.job_losers_sa L(1/`l').`var'
		estimates store m_`var'_`l'
        testparm L(1/`l').`var'
		
		// Models with Covid Mean Shift
		quietly regress D.job_losers_sa L(1/`l').D.job_losers_sa L(1/`l').`var' covid_shock
		estimates store m_`var'_`l'_cov
        testparm L(1/`l').`var' covid_shock
    }
}

/*
	Model Selection
______________________*/

estimates stats _all

** No Unit Root Assumption
arima job_losers_sa covid_shock, arima(3,0,0)
estimates store Final_Model_No_Unit_Root
estimates stats

** Unit Root Assumption // Differencing
regress D.job_losers_sa L(1/3).D.job_losers_sa L(1/3).corporate covid_shock
// Keeping only SS corp lag coefs
gen d_job_losers = D.job_losers_sa
regress d_job_losers L(1/2).d_job_losers L(1/2).corporate covid_shock
predict res, residual
pac res
estimates store Final_Model_Unit_Root
estimates stats

/*
	Forecasting
__________________*/

set seed 2004
tsappend, add(12)

replace covid_shock = 0 if missing(covid_shock)

forecast create model, replace
forecast estimates Final_Model_No_Unit_Root
forecast solve, simulate ( errors betas , statistic (stddev , prefix (sd_)) reps(1000))

gen L_int = f_job_losers_sa + invnormal(0.025)*sd_job_losers_sa
gen U_int = f_job_losers_sa + invnormal(0.975)*sd_job_losers_sa

tsline f_job_losers_sa job_losers_sa L_int U_int if t >= ym(2022, 1), legend(order(1 2 3) ///
label(1 "Point") label(2 "Actual") label(3 "95% interval") ) ///
lcolor(blue black red red) lpattern(dash solid dash dash) ///
xtitle("Time (Months)") title("Job Losers as a Percent of Total Unemployed")
graph export forecast.png, replace

capture drop f_* sd_* L_int* U_int*
// Assumption, no change
drop if _n > 704
tsappend, add(12)

replace covid_shock = 0 if missing(covid_shock)
replace corporate = L.corporate if missing(corporate)

forecast create model, replace
forecast estimates Final_Model_Unit_Root
forecast solve, simulate ( errors betas , statistic (stddev , prefix (sd_)) reps(1000))

gen L_int = f_d_job_losers + invnormal(0.025)*sd_d_job_losers
gen U_int = f_d_job_losers + invnormal(0.975)*sd_d_job_losers

tsline f_d_job_losers d_job_losers L_int U_int if t >= ym(2022, 1), legend(order(1 2 3) ///
label(1 "Point") label(2 "Actual") label(3 "95% interval") ) ///
lcolor(blue black red red) lpattern(dash solid dash dash) ///
xtitle("Time (Months)") title("Change in Job Losers as a Percent of Total Unemployed") ///
subtitle("Assumed Constant Corporate Spread")
graph export forecast_constant.png, replace


capture drop f_* sd_* L_int* U_int*
// Assumption, increasing corporate spread
drop if _n > 704
tsappend, add(12)
replace covid_shock = 0 if missing(covid_shock)
//replace corporate = L.corporate + 0.5 if missing(corporate)
forvalues i = 1/12 {
    replace corporate = L.corporate + 0.1 + rnormal(0, 0.05) if missing(corporate)
}

forecast create recession_model, replace
forecast estimates Final_Model_Unit_Root
forecast solve, simulate ( errors betas , statistic (stddev , prefix (sd_)) reps(1000))

gen L_int = f_d_job_losers + invnormal(0.025)*sd_d_job_losers
gen U_int = f_d_job_losers + invnormal(0.975)*sd_d_job_losers

tsline f_d_job_losers d_job_losers L_int U_int if t >= ym(2022, 1), legend(order(1 2 3) ///
label(1 "Point") label(2 "Actual") label(3 "95% interval") ) ///
lcolor(blue black red red) lpattern(dash solid dash dash) ///
xtitle("Time (Months)") title("Change in Job Losers as a Percent of Total Unemployed") ///
subtitle("Assumed Increasing Corporate Spread")
graph export forecast_recession.png, replace