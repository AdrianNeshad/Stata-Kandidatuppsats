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
* Spara ursprungliga data
preserve

* Skapa tidsperioder (5-årsperioder för stabilitet)
gen period = floor((year-1995)/5) + 1
label define period_lbl 1 "1995-1999" 2 "2000-2004" 3 "2005-2009" ///
    4 "2010-2014" 5 "2015-2019" 6 "2020-2024"
label values period period_lbl

* Kollapsa data till genomsnitt per period och land
collapse (mean) avg_growth = GDP_g_pc ///
         (first) initial_ln_GDP = ln_GDP_pc ///
         (mean) avg_POP_g = POP_g ///
         (mean) avg_asinh_FDIi = asinh_FDIi ///
         (mean) avg_asinh_ICP = asinh_ICP ///
         (first) COLd = COLd ///
         (first) ISO3 = ISO3, ///
         by(country_id period)

* Ta bort observationer med missing värden
drop if missing(avg_growth, initial_ln_GDP)

* ALTERNATIV 1: Grundläggande absolut konvergens (endast initial BNP)
asdoc regress avg_growth initial_ln_GDP, ///
    save(Convergence_Correct.doc) replace title(Absolut_Konvergens_Alla) ///
    label

* ALTERNATIV 2: Villkorlig konvergens med kontrollvariabler
asdoc regress avg_growth initial_ln_GDP avg_POP_g avg_asinh_FDIi avg_asinh_ICP, ///
    save(Convergence_Correct.doc) append title(Villkorlig_Konvergens_Alla) ///
    label

* ALTERNATIV 3: Konvergens för olika grupper (brittiska vs franska kolonier)
asdoc regress avg_growth initial_ln_GDP avg_POP_g avg_asinh_FDIi avg_asinh_ICP ///
    if COLd==0, ///
    save(Convergence_Correct.doc) append title(Villkorlig_Konvergens_Brittiska) ///
    label

asdoc regress avg_growth initial_ln_GDP avg_POP_g avg_asinh_FDIi avg_asinh_ICP ///
    if COLd==1, ///
    save(Convergence_Correct.doc) append title(Villkorlig_Konvergens_Franska) ///
    label

* ALTERNATIV 4: Testa om konvergenshastigheten skiljer sig mellan grupper
gen COLd_initial = COLd * initial_ln_GDP
asdoc regress avg_growth initial_ln_GDP COLd COLd_initial avg_POP_g avg_asinh_FDIi avg_asinh_ICP, ///
    save(Convergence_Correct.doc) append title(Konvergens_Gruppinteraktion) ///
    label

* Test av signifikant skillnad i konvergenshastighet
test COLd_initial = 0

* Beräkna konvergenshastighet (lambda) om koefficienten är negativ
* lambda = -ln(1 + β*T)/T, där T är antal år per period (5 år)
* Half-life = ln(2)/lambda
display "Om beta-koefficienten är negativ, beräkna konvergenshastighet:"
display "lambda = -ln(1 + beta*5)/5"
display "Half-life = ln(2)/lambda"

* Återställ ursprungliga data
restore

*** TOLKNING ***
* Negativt och signifikant β (koefficient på initial_ln_GDP) = konvergens
* Ju mer negativt β, desto snabbare konvergens
* Jämför β mellan COLd=0 och COLd=1 för att se skillnader
* COLd_initial-koefficienten testar om konvergenshastigheten skiljer sig signifikant


***Grafer
***ln_GDP_pc
preserve
collapse (mean) ln_GDP_pc, by(year COLd)
twoway (line ln_GDP_pc year if COLd == 0, lcolor(blue)) (line ln_GDP_pc year if COLd == 1, lcolor(red)), legend(label(1 "Brittisk koloni") label(2 "Fransk koloni")) title("") ytitle("Log BNP per capita") xtitle("År") name(graph_ln_GDP_pc_time, replace)
graph export "ln_GDP_pc_over_time_by_COLd.png", replace
restore

**GDP_pc
preserve
collapse (mean) GDP_pc, by(year COLd)
twoway (line GDP_pc year if COLd == 0, lcolor(blue)) (line GDP_pc year if COLd == 1, lcolor(red)), legend(label(1 "Brittisk koloni") label(2 "Fransk koloni")) title("") ytitle("BNP per capita") xtitle("År") name(graph_GDP_pc_time, replace)
graph export "GDP_pc_over_time_by_COLd.png", replace
restore

** GDP_g_pc_over_time
preserve
collapse (mean) GDP_g_pc, by(year COLd)
twoway (line GDP_g_pc year if COLd == 0, lcolor(blue)) ///
       (line GDP_g_pc year if COLd == 1, lcolor(red)), ///
       legend(label(1 "Brittisk koloni") label(2 "Fransk koloni")) ///
       title("") ytitle("BNP tillväxt per cap.") xtitle("År") ///
       name(graph_BNPb_time, replace)

graph export "GDP_g_pc_over_time_by_COLd.png", replace
restore

*** Graf för koloni och investeringar
preserve
collapse (mean) asinh_FDIi, by(year COLd)
twoway ///
  (line asinh_FDIi year if COLd == 0, lcolor(blue)) ///
  (line asinh_FDIi year if COLd == 1, lcolor(red)), ///
  legend(label(1 "Brittisk koloni") label(2 "Fransk koloni")) ///
  title("Genomsnittlig Foreign Invest. över tid") ///
  ytitle("asinh(Foreign Invest.)") xtitle("År") ///
  name(graph_FDIi_time, replace)
graph export "FDIi_by_Colony.png", replace
restore


preserve
collapse (mean) GDP_g_pc asinh_FDIi, by(year)
twoway ///
    (line GDP_g_pc year, lcolor(blue) lpattern(solid) yaxis(1)) ///
    (line asinh_FDIi year, lcolor(red) lpattern(dash) yaxis(2)), ///
    legend(label(1 "BNP-tillväxt per cap.") label(2 "asinh(FDI inflows)")) ///
    title("BNP-Tillväxt per cap. och Foreign Invest. över tid") ///
    xtitle("År") ///
    ytitle("ΔBNP  per cap.", axis(1)) ///
    ytitle("asinh(FDIi)", axis(2)) ///
    name(graph_GDP_g_pc_FDI_time, replace)
graph export "GDP_g_pc_FDI_over_time.png", replace
restore
