cd "/Users/adrianneshad/Documents/GitHub/Stata-Kandidatuppsats"

clear all
use Dataset.dta, clear
save Data.dta, replace

*** DESKRIPTIV STATISTIK ***
des

* Panelstruktur
egen country_id = group(ISO3)	
xtset country_id year

* Skapa logaritmerade och transformerade variabler
gen GDP_pc = GDP / POP
gen ln_GDP_pc = ln(GDP_pc)
gen GDP_g_pc = 100 * D.ln_GDP_pc

gen ln_POP = ln(POP)
gen asinh_ICP = asinh(ICP) 
gen asinh_FDIi = asinh(FDIi)

gen POP_g = 100 * D.ln_POP

* Lista variabler
local depvarA GDP_g_pc
local indvars COLd index c.COLd#c.index
local contvars POP_g asinh_FDIi asinh_ICP


*** DESKRIPTIV STATISTIK TILL WORD ***
asdoc tabstat `depvarA' COLd index `contvars', stat(N mean p50 sd min max) ///
    columns(statistics) save(Deskriptiv.doc), replace label

*** KORRELATIONER ***
asdoc pwcorr `depvarA' COLd index `contvars', star(0.05) columns(statistics) ///
    save(PairwiseKorrelation.doc), replace
	
*** BREUSCH-PAGAN TEST ***
regress `depvarA' `indvars' `contvars'
estat hettest
outreg2 using "BP_test_BNPa.doc", replace label ///
    addstat("Chi2", r(chi2), "p-value", r(p)) ctitle(BP_test)	

*** RAMSEY RESET / OVTEST ***
regress `depvarA' `indvars' `contvars'
estat ovtest
outreg2 using "OV_test_BNPa.doc", replace label ///
    addstat("F-statistic", r(F), "p-value", r(p)) ctitle(OV_test)

*** MULTIKOLLINEARITET ***
reg `depvarA' COLd index `contvars'
estat vif

*** FIXED EFFECTS – RANDOM EFFECTS – HAUSMAN ***
xtreg `depvarA' `indvars' `contvars', fe
est store fixedA
asdoc xtreg `depvarA' `indvars' `contvars', fe ///
    save(Hausman_test.doc) replace label nest cnames(FE)

xtreg `depvarA' `indvars' `contvars', re
est store randomA
asdoc xtreg `depvarA' `indvars' `contvars', re ///
    save(Hausman_test.doc) append label nest cnames(RE)

hausman fixedA randomA, sigmamore
asdoc hausman fixedA randomA, sigmamore ///
    save(Hausman_test.doc) append label
	

*** POOLED OLS & RANDOM EFFECTS JÄMFÖRELSE ***
reg `depvarA' `indvars' `contvars', vce(cluster country_id)
est store ols

xtreg `depvarA' `indvars' `contvars', re vce(cluster country_id)
est store random_ols

esttab ols random_ols using "Combined_OLS_RE.rtf", ///
    replace ///
    nocons ///
    label ///
    mtitle("Pooled OLS" "Random Effects") ///
    b(3) se(3) parentheses ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2 r2_a F chi2, ///
        fmt(0 3 3 2 2) ///
        labels("Observations" "R-squared" "Adj. R-squared" "F-statistic" "Chi-squared")) ///
    title("Combined OLS and RE for GDP per capita growth")
	
*** MISSING DATA CHECK ***
misstable summarize

*** TESTA INTERAKTIONSEFFEKTEN ***
gen COLd_index = COLd * index
reg `depvarA' COLd index COLd_index `contvars', vce(cluster country_id)
test COLd_index = 0



*** β‑konvergens – testar om fattigare länder växer snabbare än rikare ****

* 1) Skapa laggad log‑BNP om inte redan
gen L_ln_GDP_pc = L.ln_GDP_pc

* 2) Konvergens för alla länder (skapar filen)
asdoc regress GDP_g_pc L_ln_GDP_pc POP_g asinh_FDIi asinh_ICP, save(Convergence.doc) replace title(Alla)

* 3) Konvergens per grupp (lägger till i samma fil)
* För tidigare brittiska kolonier (COLd==0)
asdoc regress GDP_g_pc L_ln_GDP_pc POP_g asinh_FDIi asinh_ICP if COLd==0, save(Convergence.doc) append title(COLd0)

* För tidigare franska kolonier (COLd==1)
asdoc regress GDP_g_pc L_ln_GDP_pc POP_g asinh_FDIi asinh_ICP if COLd==1, save(Convergence.doc) append title(COLd1)

* 4) Tolkning: ett negativt och signifikant β‑värde (koefficient på L_ln_GDP_pc) indikerar konvergens. 
*    Jämför β för COLd=0 och COLd=1 för att se om en grupp konvergerar snabbare än den andra.
