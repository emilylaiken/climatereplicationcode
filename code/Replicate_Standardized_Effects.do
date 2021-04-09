
* THIS FILE REPLICATES THE STANDARDIZED EFFECTS SHOWN IN FIGURES 4 AND 5 OF HSIANG, BURKE, AND MIGUEL 2013
* For each study, we provide the code to calculate the standardized effect, as well as the code to generate
*	the residualized variables for the 12 studies that appear in Figure 2


clear all
clear matrix

set mem 2G
set matsize 10000

* replace the following line with the directory where the replication data was unzipped
global wd "\Documents\Dropbox\Marshall-Sol\drafts\Science_review_v4\replication"
//global wd "~\Dropbox\SHARED_FOLDERS\Marshall-Sol\drafts\Science_review_v4\replication" // sol's machine

* create log file
capture log close
log using $wd\output\Replicate_Standardized_Effects.txt, replace text


* work out of the code\stata_subfunctions directory because subfunctions are used in estimation
cd "$wd\code\stata_subfunctions"


* directories for datasets
global bl = "$wd\data\Bergholt_Lujala_2012\"
global bs = "$wd\data\Bohlken_Sergenti_2010\"
global bc = "$wd\data\Bruckner_Ciccone_2011\"
global burke = "$wd\data\Burke_2012\"
global burket = "$wd\data\Burke_et_al_2009\"
global bule = "$wd\data\Burke_Leigh_2010\"
global cd = "$wd\data\Card_Dahl_2011\"
global djo = "$wd\data\Dell_Jones_Olken_2012\" 
global hs = "$wd\data\Hendrix_Salehyan_2012\"
global hmc = "$wd\data\Hsiang_et_al_2011\"
global ol = "$wd\data\OLoughlin_et_al_2012\"
global mi = "$wd\data\Miguel_2005\"
global mss = "$wd\data\Miguel_Satyanath_Sergenti_2004\"
global th = "$wd\data\Theisen_2012\"
global hnn = "$wd\data\Hidalgo_Naidu_Nichter_2010\"
global jlm = "$wd\data\Jacob_Lefgren_Moretti_2007\"
global bht = "$wd\data\Buhaug_Theisen_Holtermann_2011\"
global levy = "$wd\data\Levy_et_al_2005\"
global lt = "$wd\data\Larrick_Timmerman_2011\"
global ra = "$wd\data\Ranson_2012\"
global fvu = "$wd\data\Fjelde_vonUexkull_2012\"


***** Bergholt and Lujala 2012 [76]

	use "$bl\B%26L_2012_replication_data.dta"
	xtreg onset2 climdis lpoppwt1 polity21 polity2s1 peaceyrs _spline1 _spline2 _spline3 t37-t63 if year > 1979, fe robust cluster(ccode)
	local eff = _b[climdis]
	loneway climdis ccode 
	local sd = r(sd_w)
	summ onset2 
	di `eff'/r(mean)*`sd'*100 

clear


***** Bohlken and Sergenti 2010 [44]

	use "$bs\bs_riots.dta", clear
	tsset STATE_NB YEAR 

	* re-code rain to be "rainfall loss", to match signs
	gen rain_loss = - RAIN
	gen rain_loss1 = -RAIN1

	* Run regressions and calculate standardized effects
	areg CNTVIOL rain_loss rain_loss1 i.YEAR, a(STATE_NB) cl(STATE_NB)
	local eff = _b[rain_loss]
	loneway rain_loss STATE_NB
	local sd = r(sd_w)
	summ CNTVIOL
	di `eff'/r(mean)*`sd'*100 

	* For Fig 2
	//so units in % of mean:
	sum CNTVIOL
	loc mean_riots = r(mean)
	replace CNTVIOL = CNTVIOL/`mean_riots'*100 

	// dropping extreme outlier years
	drop if STATE == "Gujarat" & YEAR == 1988
	drop if STATE == "Haryana" & YEAR == 1988

	qui reg CNTVIOL    I* RAIN1
	predict riots_r, resid
	qui reg rain_loss    I* RAIN1
	predict rain_r, resid
	outsheet riots_r rain_r using "$wd\output\DataForFigure2\BOLKEN.csv", comma replace

clear

***** Bruckner and Ciccone 2011 [79]

	use "$bc\Bruckner_Ciccone_2011.dta"

	* For standardized effects plot
	areg polity_change lgpcp_l lgpcp_l2 ccode#c.year i.year  , a(ccode) cluster(ccode)
	local eff = _b[lgpcp_l2]
	loneway lgpcp_l2 country
	local sd = r(sd_w)
	summ polity_change
	di `eff'/r(mean)*`sd'*100

clear


***** Burke 2012 [72]

	* Standardized effects
	use "$burke\Burke_2012.dta"

	xi i.ccodewb*year, pref(_yi_)
	xi i.ccodewb, pref(_i_)

	//since pre-divided by 100 in paper
	replace ins_tem1 = ins_tem1*100
	replace ins_pre = ins_pre*100

	// baseline regression

	xtreg dummy L(0/1).ins_tem1 L(0/1).ins_pre  _yi_* if sample_act==1, fe cluster(ccodewb)
	local eff = _b[ins_tem1]
	local se = _se[ins_tem1]
	loneway ins_tem1 ccodewb if sample_act==1
	local sd = r(sd_w)
	sum dummy if e(sample)
	di `eff'/r(mean)*`sd'*100 

	* Generate residualized variables and export dataset for Fig 2
	qui reg dummy L(0).ins_tem1 L(0).ins_pre  _yi_* if sample_act==1, cluster(ccodewb)
	sum dummy if e(sample)
	loc mean_dummy = r(mean)
	qui reg L(0/1).ins_tem1 L(0/1).ins_pre _yi_* if sample_act==1 & dummy ~=.
	predict resid_tem1 if e(sample), resid
	qui reg dummy L(1).ins_tem1 L(0/1).ins_pre _yi_* if sample_act==1 & ins_tem1 ~=.
	predict resid_dummy if e(sample), resid
	replace resid_dummy = resid_dummy*100/`mean_dummy' // so units are in % of mean

	//Oman and UAE (1977) have extreme high temp
	//Finland (1985, 1987) & Swedend 1985 have extreme low temp
	keep if resid_tem1 < 2 & resid_tem1 > -2
	outsheet resid_dummy resid_tem1 using "$wd\output\DataForFigure2\BURKE_LEADER.csv", comma replace

clear


***** Burke et al 2009 [65]

	use "$burket\Burke_etal_2009.dta", clear
	* Regression we report:  Model 2 in Table 1
	xtreg war_prio_new temp_all temp_all_lag  prec_all prec_all_lag Iccyear*, fe i(ccode) cl(ccode)
	* standardized effect
	qui xtreg war_prio_new temp_all temp_all_lag  prec_all prec_all_lag Iccyear*, fe i(ccode) cl(ccode)
	local eff = _b[temp_all]
	local se = _se[temp_all]
	loneway temp_all ccode if smpl<.
	local sd = r(sd_w)
	summ war_prio_new if smpl<.
	di `eff'/r(mean)*`sd'*100	

	* Residualized variables for fig 2
	qui xtreg war_prio_new temp_all_lag  prec_all prec_all_lag Iccyear*, fe i(ccode) cl(ccode)
	predict war_partial if e(sample), e
	qui xtreg temp_all temp_all_lag  prec_all prec_all_lag Iccyear*, fe i(ccode) cl(ccode)
	predict temp_partial if e(sample), e
	summ war_prio_new
	gen war_partial_pct = war_partial/`r(mean)'*100

	outsheet war_partial_pct temp_partial using "$wd\output\DataForFigure2\BURKE_PNAS.csv", comma replace

clear


***** Burke and Leigh 2010 [78]

	use "$bule\Burke_Leigh_2010.dta"
	
	areg demchangeevent L(1/2).precipitationinteract L(1/2).tempdevnew1interact  i.year if sample2001==1 & code~="LBY", absorb(ccode) cluster(ccode)
	matrix b = e(b)
	scalar eff=b[1,3]

	*statistics for std effects
	loneway tempdevnew1interact ccode if e(sample)
	local sd = r(sd_w)
	sum demchangeevent if e(sample) 
	di eff/r(mean)*`sd'*100

clear



***** Hsiang Meng Cane 2011
	use "$hmc\Cane_Hsiang_Meng_2011.dta"

	drop if year < 1950
	drop if year > 2004
	drop if year == 1989
	drop if region_dummy == 0

	* main regression
	reg conflicts_per_country_m3_pop nino3_late post1989 year if enso_region=="TE"
	local eff = _b[nino3_late]
	summ nino3_late
	local sd = r(sd_w)
	summ conflicts_per_country_m3_pop 
	di `eff'/r(mean)*`sd'*100

	* FOR PLOT
	//so units in % of mean:
	sum conflicts_per_country_m3_pop
	loc mean_ACR_TE = r(mean)
	replace conflicts_per_country_m3_pop = conflicts_per_country_m3_pop/`mean_ACR_TE'*100 

	//so units in % of mean:
	sum conflicts_per_country_global
	loc mean_ACR_ALL = r(mean)
	replace conflicts_per_country_global= conflicts_per_country_global/`mean_ACR_ALL'*100

	qui reg conflicts_per_country_m3_pop  post1989 year
	predict conflicts_r, resid
	replace conflicts_r = conflicts_r

	qui reg conflicts_per_country_global  post1989 year
	predict conflicts_global_r, resid
	replace conflicts_global_r = conflicts_global_r

	reg nino3_late  post1989 year
	predict nino_r, resid
	outsheet conflicts_r nino_r using "$wd\output\DataForFigure2\HSIANG.csv", comma replace

clear	


***** Card and Dahl 2011 [37]

	* Prepare data, following their merge script

	clear all
	*change just to work with Card and Dahl data
	cd $cd 

	***Merge nfl and nibrs data together
	*Note: Must first create nibrsihour dataset using their crimedata-step1.do, crimedata-step2.do, and crimedata-step3.do		
	use nibrsihour, replace
	merge teama mdy using nfl-online
	replace dow=dow(mdy)
	drop _merge
	sort statefips mdy
	save nibrsnfl, replace

	***Merge in holidays & weather for online appendix
	merge statefips mdy using holidaysweather
	tab _merge
	drop if _merge!=3
	drop _merge
	compress

	***create cumweek (week of season) for non-gameday Sundays in our sample
	*Note, due to Christmas falling on a Sunday in 2005, no Sunday games (for the teams in our sample, although a few NFL games were played) in the 16th week of the season
	replace cumweek=16 if mdy==mdy(12,25,2005) 
	keep if cumweek!=.
	collapse cumweek, by(mdy) fast
	rename cumweek newcumweek
	sort mdy
	save mdytocumweek, replace
	use all, replace
	sort mdy
	merge mdy using mdytocumweek
	drop _merge
	drop cumweek
	rename newcumweek cumweek
	gen season=year
	replace season=year-1 if month==1|month==2

	*create season (yy) and week of season (ww) indicators
	tab season, gen(yy)
	tab cumweek, gen(ww)

	*Create variables based on pre-game spread and halftime spread
	gen spreadcat=1*(spread<-7 & spread!=.) + 1*(spread>=-7 & spread<-3.5) + 2*(spread>=-3.5 & spread<=3.5) + 3*(spread>3.5 & spread<=7) + 3*(spread>7 & spread!=.)
	gen byte upsetloss=(gameday & loss & spreadcat==1)
	gen byte closeloss=(gameday & loss & spreadcat==2)
	gen byte upsetwin=(gameday & win & spreadcat==3)
	gen byte predwin=spreadcat==1
	gen byte predclose=spreadcat==2
	gen byte predloss=spreadcat==3
	gen spreadmiss=spread==.
	replace spread=0 if spread==.
	gen halfspread=-halfdiff
	gen halfpredwin=halfspread<-3
	gen halfpredclose=halfspread>=-3 & halfspread<=3
	gen halfpredloss=halfspread>3 & halfspread!=.
	gen halfupsetloss=halfpredwin & loss
	gen halfcloseloss=halfpredclose & loss
	gen halfupsetwin=halfpredloss & win

	replace halfspread=0 if halfspread==.

	gen orcombo=(riv | s4t4p80) & !lowplayoff

	global cat "riv lowplayoff s4t4p80 orcombo"

	foreach var of varlist $cat {
	  gen byte upsetloss`var'=upsetloss*`var'==1
	  gen byte closeloss`var'=closeloss*`var'==1
	  gen byte upsetwin`var'=upsetwin*`var'==1
	  gen byte predwin`var'=predwin*`var'==1
	  gen byte predclose`var'=predclose*`var'==1
	  gen byte predloss`var'=predloss*`var'==1
	}

	gen byte upsetlossbase=upsetloss
	gen byte closelossbase=closeloss
	gen byte upsetwinbase=upsetwin
	gen byte predwinbase=predwin
	gen byte predclosebase=predclose
	gen byte predlossbase=predloss

	foreach i in 1 4 {
	  gen upsetloss`i'=upsetloss & easterntime==`i'
	  gen closeloss`i'=closeloss & easterntime==`i'
	  gen upsetwin`i'=upsetwin & easterntime==`i'
	  gen predwin`i'=predwin & easterntime==`i'
	  gen predclose`i'=predclose & easterntime==`i'
	  gen predloss`i'=predloss & easterntime==`i'
	}

	* THIS NEXT CODE IS FROM THEIR tables.do SCRIPT

	*create a sample for estimating the robustness check on how to treat missings as zeros
	gen insamp=regularseason & dow==0 & _fillintot==0

	*create an estimation sample which only includes observations in baseline regression
	global basevars="upsetlossbase closelossbase upsetwinbase predwinbase predclosebase predlossbase"
	global smallhh = "christeve christday newyeareve newyearday halloween thankswkd laborwkd columwkd vetwkd"
	global weather="hot hiheatindx cold windy anyrain anysnow"
	xtpoisson ipmfhometot $basevars ww* yy* $smallhh $weather if insamp, fe i(orinum)
	keep if e(sample)
	egen teamseason = group(teama season)
	qui tab orinum, gen(oo)

	* FINALLY, RUN, THE REGRESSIONS
	* with climate variables only. this is very slow.
	poisson ipmfhometot hot cold ww* yy* oo* if insamp, difficult iterate(25) cluster(teamseason)  
	local eff = _b[hot]
	loneway hot orinum 
	local sd = r(sd_w)
	di `eff'*`sd'*100

	*without clustering. much faster because we can use xtpoisson
	xtpoisson ipmfhometot hot cold ww* yy* if insamp, fe i(orinum)  //no clustering

clear	

*return to working directory
cd "$wd\code\stata_subfunctions" 
	
	

***** Dell Jones Olken 2012 [21] (both conflict onset and irregular leader transition)

	*first construct variables following their code and then calculate standardized effects at bottom
	use "$djo\climate_panel.dta", clear
	global rfe = 1 //1 for region*year, 2 for year only
	global maineffectsonly = 0 //1 to drop all interactions

	keep if year <= 2003
	encode parent, g(parent_num)
	gen lgdp=ln(gdpLCU)
	encode fips60_06, g(cc_num)
	sort country_code year
	tsset cc_num year

	* calculate GDP growth
	gen temp1 = l.lgdp
	gen g=lgdp-temp1
	replace g = g * 100 
	drop temp1
	summarize g

	g lnag = ln(gdpWDIGDPAGR) 
	g lnind = ln(gdpWDIGDPIND) 
	g lninvest = ln(rgdpl*ki/100)
	g lngdpwdi = ln(gdpLCU)
	foreach X in ag ind gdpwdi invest {
		g g`X' = (ln`X' - l.ln`X')*100
	}

	* Drop if less than 20 yrs of GDP data
	g tempnonmis = 1 if g != .
	replace tempnonmis = 0 if g == .
	bys fips60_06: egen tempsumnonmis = sum(tempnonmis)
	drop if tempsumnonmis  < 20

	* Make sure all subcomponents are non-missing in a given year
	g misdum = 0
	for any ag ind : replace misdum = 1 if gX == .
	for any ag ind : replace gX = . if misdum == 1

	preserve
	keep if lnrgdpl_t0 < . 
	bys fips60_06: keep if _n == 1 
	xtile initgdpbin = ln(lnrgdpl_t0), nq(2)
	keep fips60_06 initgdpbin
	tempfile tempxtile
	save `tempxtile',replace
	restore

	merge n:1 fips60_06 using `tempxtile'
	drop _merge
	tab initgdpbin, g(initxtilegdp)

	preserve
	keep if wtem50 < . 
	bys fips60_06: keep if _n == 1 
	xtile initwtem50bin = wtem50 , nq(2)
	keep fips60_06 initwtem50bin
	save `tempxtile',replace
	restore

	merge n:1 fips60_06 using `tempxtile'
	drop _merge
	tab initwtem50bin, g(initxtilewtem)
	preserve
	keep if year == 1995
	sort fips60_06 year
	by fips60_06: keep if _n == 1
	g temp = gdpSHAREAG 
	xtile initagshare1995 = ln(temp), nq(2)
	replace initagshare1995 = . if gdpSHAREAG == .
	keep fips60_06 initagshare1995 
	tempfile tempxtile
	save `tempxtile',replace
	restore
	
	merge n:1 fips60_06 using `tempxtile'
	drop _merge
	tab initagshare1995 , g(initxtileagshare)

	tsset	
	foreach Y in wtem wpre  {
		gen `Y'Xlnrgdpl_t0 =`Y'*lnrgdpl_t0 
		for var initxtile*: gen `Y'_X =`Y'*X
		label var `Y'Xlnrgdpl_t0 "`Y'.*inital GDP pc"
		for var initxtile*: label var `Y'_X "`Y'* X"
	}

		
	* make lags
	capture {
		for var wtem* wpre*: g fdX = X - l.X \ label var fdX "Change in X"
		for var wtem* wpre*: g L1X = l1.X 
		for var wtem* wpre*: g L2X = l2.X 
		for var wtem* wpre*: g L3X = l3.X 
		for var wtem* wpre*: g L4X = l4.X 
		for var wtem* wpre*: g L5X = l5.X 
		for var wtem* wpre*: g L6X = l6.X 
		for var wtem* wpre*: g L7X = l7.X 
		for var wtem* wpre*: g L8X = l8.X 
		for var wtem* wpre*: g L9X = l9.X 
		for var wtem* wpre*: g L10X = l10.X
	}	

	tab year, gen (yr)
	local numyears = r(r) - 1
	if $rfe == 1 {
		foreach X of num 1/`numyears' {
				foreach Y in MENA SSAF LAC WEOFF EECA SEAS {
					quietly gen RY`X'X`Y'=yr`X'*_`Y'
					quietly tab RY`X'X`Y'
				}
				quietly gen RYPX`X'=yr`X'*initxtilegdp1
			}
	}
	else if $rfe == 2 {
		foreach X of num 1/`numyears' {
				quietly gen RY`X'=yr`X'
			}
	}

	
	* Make the political variables, and then merge into the main dataset
	tempfile initialdata temppolity tempfips temparchigos tempconflict 
	* drop _merge
	sort fips60_06 year
	save `initialdata', replace
	use "$djo\polity", clear
	save `temppolity', replace

	* Code Archigos data
	clear
	insheet using "$djo\fips_to_cow.csv"
	sort ccode
	save `tempfips', replace
	use "$djo\Archigos.dta", clear

	*--create some variables--*
	replace startdate="09/11/1942" if startdate=="9/11/1942"
	replace startdate="01/01/2001" if startdate=="1/1/2001"
	replace startdate="01/01/2002" if startdate=="1/1/2002"
	replace startdate="14/09/1923" if startdate=="14/9/1923"
	keep ccode startdate entry
	gen year=substr(startdate,7,4)
	destring year, replace
	drop if year<1950
	tab year
	recode entry 2=1 //transition imposed by other country coded as irregular
	tab entry
	sort ccode

	replace ccode=679 if ccode==680 //combining Yemen pre and post 1990
	replace ccode=255 if (ccode==260 & year==1998) //separating post-unification Germany from West Germany
	sort ccode
	merge ccode using `tempfips', uniqusing
	tab _merge
	keep if _merge==3 //_merge==1s are places like Bavaria that no longer exist and East/West Germany, South Vietnam, which we are dropping and thus were excluded from the code conversion list without leader data
		//merge==2s are places in the code list without transition data, like Monaco, Kiribat, and the Federal States of Micronesia
	drop _merge

	*---code countries whose borders have changed---*
	*aggregate to country*year observations
	rename fips_code fips
	gen fips60_06=fips
	replace fips60_06=	"PBD"	if (fips==	"PK"	& year <	1971	)
	replace fips60_06=	"ETR"	if (fips==	"ET"	& year <	1993	)
	replace fips60_06=	"NSA"	if (fips==	"SF"	& year <	1990	)
	replace fips60_06=	"SVT"	if (fips==	"RS"	& year <	1992	)

	//codebook takes care of Yugoslavia and Czechoslovakia, because these have COW codes
	sort fips60_06

	*---create some variables---*
	collapse (sum)entry, by (year fips60_06)
	recode entry 2/7=1 //some countries have more than one transition in a given year. If any of the transitions were irregular, I code that year as experiencing an irregular transition
	sort fips60_06 
	save `temparchigos', replace

	*Conflict code
	clear
    insheet using "$djo\fips_to_cow.csv"
    sort ccode
    save temp, replace

	clear
	insheet using "$djo\prio.csv"
	rename gwnoa countrya 
	rename gwnob countryb
	gen countryc=.
	gen countryd=.
	gen countrye=.
	gen countryf=.
	replace countryc=666 if countrya=="220, 666, 200"
	replace countryc=645 if countryb=="651, 645, 663, 660, 652"
	replace countryc=200 if countryb=="900, 200, 2"
	replace countryd=200 if countrya=="220, 666, 200"
	replace countryd=663 if countryb=="651, 645, 663, 660, 652"
	replace countryd=2 if countryb=="900, 200, 2"
	replace countrye=660 if countryb=="651, 645, 663, 660, 652"
	replace countryf=652 if countryb=="651, 645, 663, 660, 652"
	replace countrya="220" if countrya=="220, 666, 200"
	replace countryb="651" if countryb=="651, 645, 663, 660, 652"
	replace countryb="900" if countryb=="900, 200, 2"
	destring countrya, replace
	destring countryb, replace
	replace countrya=. if countrya==-99
	g intensitycivil = intensity

	replace intensitycivil = 0 if type != 3 & type != 4 & intensitycivil != .

	gen obs=_n
	reshape long country, i(obs) j(temp) str
	keep year intensity country intensitycivil
	drop if country==.
	drop if year<1950

	rename country ccode
	replace ccode=679 if ccode==680 //combining Yemen pre and post 1990
	sort ccode
	merge ccode using temp, uniqusing
	tab _merge
	//the _merge==1s are Yemen (pre-1950, before climate data) and South Vietnam
	//the _merge==2s are countries not in the conflict table
	keep if _merge==3
	drop _merge

	*---code countries whose borders have changed---*

	*aggregate to country*year observations
	rename fips_code fips
	gen fips60_06=fips
	replace fips60_06=	"PBD"	if (fips==	"PK"	& year <	1971	)
	replace fips60_06=	"ETR"	if (fips==	"ET"	& year <	1993	)
	replace fips60_06=	"NSA"	if (fips==	"SF"	& year <	1990	)
	replace fips60_06=	"SVT"	if (fips==	"RS"	& year <	1992	)

	//codebook takes care of Yugoslavia and Czechoslovakia, because these have COW codes

	* Merge conflict and create vars
	collapse intensity intensitycivil, by (fips60_06 year)
    sort fips60_06 year
    merge fips60_06 year using `initialdata', unique
	tab _merge

	//the _merge==1s are Georgia, Croatia, and Bosnia. They are in the conflicts data before I code them as a state (i.e. conflict in Yugoslavia -> split categorized as conflict in Croatia and Bosnia)
	drop if _merge==1
	drop _merge
	replace intensity=0 if intensity==.
	tab intensity

	collapse (max) intensity, by (fips60_06 year)
	sort fips60_06 year
    merge fips60_06 year using `initialdata', unique
	tab _merge

	//the _merge==1s are Georgia, Croatia, and Bosnia. They are in the conflicts data before I code them as a state (i.e. conflict in Yugoslavia -> split categorized as conflict in Croatia and Bosnia)
	drop if _merge==1
	drop _merge
	sort cc_num year
	tsset cc_num year

	gen lintensity=l.intensity
	drop if year==1950 
	gen change=.
	replace change=1 if lintensity!=intensity
	recode change .=0
	gen mod2n=.
	replace mod2n=1 if (intensity==0 & lintensity==1)
	replace mod2n=0 if (intensity!=0 & lintensity==1)
	replace mod2n=. if (lintensity==0 | lintensity==2)
	gen int2n=.
	replace int2n=1 if (intensity==0 & lintensity==2)
	replace int2n=0 if (intensity!=0 & lintensity==2)
	replace int2n=. if (lintensity==0 | lintensity==1)
	gen war2n=.
	replace war2n=1 if (intensity==0 & (lintensity==1 | lintensity==2))
	replace war2n=0 if (intensity!=0 & (lintensity==1 | lintensity==2))
	replace war2n=. if lintensity==0
	gen n2mod=.
	replace n2mod=1 if (intensity==1 & lintensity==0)
	replace n2mod=0 if (intensity!=1 & lintensity==0)
	replace n2mod=. if lintensity!=0
	gen n2int=.
	replace n2int=1 if (intensity==2 & lintensity==0)
	replace n2int=0 if (intensity!=2 & lintensity==0)
	replace n2int=. if lintensity!=0
	gen n2war=.
	replace n2war=1 if (intensity==1 & lintensity==0)
	replace n2war=1 if (intensity==2 & lintensity==0)
	replace n2war=0 if (intensity==0 & lintensity==0)
	replace n2war=. if lintensity!=0



	*-----------------------gen a parent country var for clustering------------------------------------*
	capture{
		gen parentconflict=fips
		replace parentconflict=	"PBD"	if (fips==	"PK"	)
		replace parentconflict=	"PBD"	if (fips==	"BG"	)
		replace parentconflict=	"CZK"	if (fips==	"LO"	)
		replace parentconflict=	"CZK"	if (fips==	"EZ"	)
		replace parentconflict=	"ETR"	if (fips==	"ET"	)
		replace parentconflict=	"ETR"	if (fips==	"ER"	)
		replace parentconflict=	"YGL"	if (fips==	"BK"	)
		replace parentconflict=	"YGL"	if (fips==	"HR"	)
		replace parentconflict=	"YGL"	if (fips==	"MK"	)
		replace parentconflict=	"YGL"	if (fips==	"SI"	)
		replace parentconflict=	"YGL"	if (fips==	"YI"	)
		replace parentconflict=	"NSA"	if (fips==	"SF"	)
		replace parentconflict=	"NSA"	if (fips==	"WA"	)
		replace parentconflict=	"SVT"	if (fips==	"AM"	)
		replace parentconflict=	"SVT"	if (fips==	"AJ"	)
		replace parentconflict=	"SVT"	if (fips==	"BO"	)
		replace parentconflict=	"SVT"	if (fips==	"EN"	)
		replace parentconflict=	"SVT"	if (fips==	"GG"	)
		replace parentconflict=	"SVT"	if (fips==	"KZ"	)
		replace parentconflict=	"SVT"	if (fips==	"KG"	)
		replace parentconflict=	"SVT"	if (fips==	"LG"	)
		replace parentconflict=	"SVT"	if (fips==	"LH"	)
		replace parentconflict=	"SVT"	if (fips==	"MD"	)
		replace parentconflict=	"SVT"	if (fips==	"RS"	)
		replace parentconflict=	"SVT"	if (fips==	"TI"	)
		replace parentconflict=	"SVT"	if (fips==	"TX"	)
		replace parentconflict=	"SVT"	if (fips==	"UP"	)
		replace parentconflict=	"SVT"	if (fips==	"UZ"	)
		}
		merge 1:1 fips60_06 year using `temparchigos'

	tab _merge
	drop if _merge==2 
	gen lt=.
	replace lt=0 if _merge==1
	replace lt=1 if _merge==3
	tab lt
	count
	gen regtr=.
	replace regtr=1 if entry==0
	recode regtr .=0
	gen irregtr=.
	replace irregtr=1 if entry==1
	recode irregtr .=0
	drop _merge

	*--Create a region x year variable for clustering
	g region=""
	foreach X in _MENA   _SSAF   _LAC    _WEOFF  _EECA   _SEAS {
		replace region="`X'" if `X'==1
	}
	g regionyear=region+string(year)
	encode regionyear, g(rynum)	

	* FINALLY, RUN REGRESSIONS WE WANT AND CALCULATE STANDARDIZED EFFECTS
	gen poor=(initgdpbin==1)  //indicator for poor

		* Irregular Leader Transition
		cgmreg irregtr wtem wtem_initxtilegdp1 wpre wpre_initxtilegdp1  RY* i.cc_num, cluster(parent_num rynum)
		lincom wtem + wtem_initxtilegdp1  //this is the effect in poor countries
		local eff = r(estimate)
		local se = r(se)
		loneway wtem fips if poor
		local sd = r(sd_w)
		summ irregtr if poor
		di `eff'/r(mean)*`sd'*100

		* Conflict Onset
		cgmreg n2war wtem wtem_initxtilegdp1 wpre wpre_initxtilegdp1  yr* i.cc_num, cluster(parent_num rynum)
		lincom wtem + wtem_initxtilegdp1  //this is the effect in poor countries
		local eff = r(estimate)
		local se = r(se)
		loneway wtem fips if poor
		local sd = r(sd_w)
		summ n2war if poor
		di `eff'/r(mean)*`sd'*100

clear


***** Hendrix and Selahyan 2012 [46]

	use "$hs\Hendrix_Salehyan_JPR491_Replication.dta"
	gen l_violent_events = log(violent_events)
	gen abs_precip_dev = abs(GPCP_precip_mm_deviation_sd)

	* regression we report
	areg l_violent_events abs_precip_dev i.year, a(ccode) cl(ccode) 
	local eff = _b[abs_precip_dev]

	* For standardized effects. precip is already standardized
	summ l_violent_events
	di `eff'/r(mean)

	* For Fig 2 plot.
	xi i.ccode, pref(_i_)
	xi i.year, pref(_y_)

	qui reg l_violent_events  _i_* _y_* if abs_precip_dev ~ = .
	predict resid_l_violent_events if e(sample), resid
	qui reg abs_precip_dev  _i_* _y_* if l_violent_events ~= .
	predict resid_abs_precip_dev if e(sample), resid
	replace resid_l_violent_events = resid_l_violent_events*100
	drop if resid_l_violent_events == . & resid_abs_precip_dev == .
	outsheet resid_l_violent_events resid_abs_precip_dev using "$wd\output\DataForFigure2\HENDRIX.csv", comma replace

clear



***** Miguel 2005 [40]

	use "$mi\Miguel_2005.dta"
	xi i.vid2*year, pref(_yi_)
	xi i.vid2, pref(_i_)
	xi i.year, pref(_y_)

	* main regression we report: column 5 of table 4
	areg witch_murders any_rain _y_* [aw=kaya], a(vid2) cluster(vid2)
	local eff = _b[any_rain]
	loneway any_rain vid2 [aw=kaya] if e(sample)
	local sd = r(sd_w)
	summ witch_murders [aw=kaya] if e(sample)
	di `eff'/r(mean)*`sd'*100

clear

	

***** Miguel Satyanath Sergenti 2004 

	use "$mss\mss_repdata.dta"

	* main regression: country FE and country time trends
	areg any_prio GPCP_g GPCP_g_l ccode#c.year, a(ccode) cluster(ccode)

	* standardized effects
	local eff = _b[GPCP_g_l]
	local se = _se[GPCP_g_l]
	loneway GPCP_g ccode
	local sd = r(sd_w)
	summ any_prio
	di `eff'/r(mean)*`sd'*100

clear



***** O'Loughlin et al 2012

	insheet using "$ol\longDataMonth.csv"
	gen year = floor(mo_yr/100)
	gen mo = mo_yr-year*100
	egen grid_mo = group(grid_id mo)

	* For standardized effect
	* Twoway clustering at grid and country-year
	egen cty_yr = group(country year)
	xi i.year  //have to do it like this bc xtivreg does not allow factor variables; using xtivreg2 because it allows multi-way clustering
	xtivreg2 events spi6 ti6 _I*, fe i(grid_mo) cluster(grid_id cty_yr) small

	* calculate standardized effect
	local eff = _b[ti6]
	loneway ti6 grid_mo
	local sd = r(sd_w)
	summ events
	loc mean_conflict = r(mean)
	di `eff'/`mean_conflict'*`sd'*100

	* output for Figure 2
	qui areg events spi6 i.year, a(grid_mo) 
	predict yresid, resid
	replace yresid = yresid * 100 / `mean_conflict' // so units are in % of mean
	
	qui areg ti6 spi6 i.year, a(grid_mo) 
	predict t_resid, resid
	reg yresid t_resid
	outsheet events yresid t_resid grid_id year mo cty_yr using "$wd\output\DataForFigure2\OLAUGHLIN.csv", comma replace
	

clear


***** Theisen 2012 [14]

	// To run this code, user needs subfunction ols_spatial_HAC.ado, described and available at http://www.solomonhsiang.com/computing/stata-code
	//	This .ado calculates Conley standard errors

	use "$th\Theisen_2012.dta"

	* baseline regression, running with spatial standard errors
	gen lon = column/4  //his grids are 0.25deg on a side, so divide by 4 to get in units of degrees
	gen lat = row/4

	//to run spatial code, need to manually construct all variables, including constant
	gen constant = 1
	xi i.et_id, pref(_i_)
	xi i.year, pref(_y_)

	// an infinite lag length is equal to clustering at the grid level
	// "dropvar" is needed to drop collinear vars
	// "bartlett" uses a linear bartlett kernel in space (a cone), otherwise it is uniform.

	ols_spatial_HAC conflict temp prec constant _y_* _i_*, lat(lat) lon(lon) p(et_id) t(year) lag(1000) dist(100) dropvar bartlett
	reg conflict temp prec _y_* _i_*
	local eff = _b[temp]
	local se = _se[temp]
	loneway temp et_id if e(sample)
	local sd = r(sd_w)
	sum conflict if e(sample)
	di `eff'/r(mean)*`sd'*100

		
	* For Fig 2 plot
	sum conflict if e(sample)
	loc mean_conflict = r(mean)
	qui reg conflict prec  y1990  y1991 y1992 y1993 y1994 y1995 y1996 y1997 y1998 y1999 y2000 y2001 y2002 y2003 y2004 _i_* if temp ~= .
	predict resid_conflict_onset if e(sample), resid
	qui reg  temp prec  y1990  y1991 y1992 y1993 y1994 y1995 y1996 y1997 y1998 y1999 y2000 y2001 y2002 y2003 y2004 _i_* if conflict ~= .
	predict resid_temp if e(sample), resid
	replace resid_conflict_onset = resid_conflict_onset * 100 /`mean_conflict' // so units are in % of mean
	outsheet resid_conflict_onset resid_temp using "$wd\output\DataForFigure2\THEISEN.csv", comma replace

clear


***** Hidalgo Naidu Nichter 2010

	use "$hnn\Hidalgo_Naidu_Nichter_2010.dta"

	* regression we report
	areg occs rain_yearly i.year, a(code) cluster(code)

	* stuff for standardized effects plot
	local eff = _b[rain_yearly]
	local se = _se[rain_yearly]
	loneway rain_yearly code if occs~=.
	local sd = r(sd_w)
	sum occs if rain_yearly ~= .
	di `eff'/r(mean)*`sd'*100

	* For Fig2 plot
	gen occupy_prob = occs
	replace occupy_prob = 1 if occs > 0 & occs ~= .
	sum occs if rain_yearly ~= .
	loc mean_occs = r(mean)
	replace occs = occs/`mean_occs'*100 if occs ~= . //so units are in % of mean
	
	xi i.year, pref(_y_)
	qui areg occs  _y_* if rain_yearly ~= ., absorb(code) cluster(code)
	predict occs_predict if e(sample), xbd
	gen occs_resid = occs - occs_predict
	qui areg  rain_yearly _y_* if occs ~= ., absorb(code) cluster(code)
	predict rain_predict if e(sample), xbd
	gen rain_resid = rain_yearly - rain_predict

	outsheet occs_resid rain_resid using "$wd\output\DataForFigure2\HIDALGO.csv", comma replace


clear


***** Jacob Lefgren Moretti 2007

	use "$jlm\Jacob_Lefgren_Moretti_2007.dta"

	* FEs:  month, jurisdiction x year, and jurisdiction specific fourth order polynomial in d.o.y
	* Construct data following their code

	*gen polynomial in d.o.y.
	forvalues i = 2/4 {
		gen doy`i'= doy^`i'
	}

	*indicator for jurisdictions
	egen juris = group(city)  //jurisdictions, presumably

	*state-by-year indicator
	egen state_yr = group(state_fips year)

	* normalize the violence and prop crime measures.
	# delimit;
	for var viol prop:
		gen ln_X=ln(X) \
		label variable ln_X "log crimes" \
		egen avg_X=mean(X), by(ori) \
		label variable avg_X "average 7 day incidence of crime" \
		gen norm_X=X/avg_X \
		label variable norm_X "crimes normalized by average frequency";
	#delimit cr
	
	* generate weights
	gen avg_tot = avg_viol + avg_prop

	* Violent crime
	cgmreg norm_viol tmean prcp i.month juris#(c.doy##c.doy##c.doy##c.doy) juris#i.year [aw=avg_tot], cl(juris state_yr)
	summ norm_viol [aw=avg_tot]  //should be equal to ~1, because it's normalized

	* get SD of temperature variable for standardized effects
	egen jurisyr = group(juris year)
	xtset jurisyr
	qui xtreg tmean juris#c.doy juris#c.doy2 juris#c.doy3 juris#c.doy4 [aw=avg_tot], fe
	predict tresid, e
	summ tresid

	* Property Crime
	cgmreg norm_prop tmean prcp i.month juris#(c.doy##c.doy##c.doy##c.doy) juris#i.year [aw=avg_tot], cl(juris state_yr)
	summ norm_prop [aw=avg_tot]

	*For Fig 2 plot
	xtset jurisyr
	qui xtreg norm_viol prcp i.month juris#c.doy juris#c.doy2 juris#c.doy3 juris#c.doy4  if tmean ~= . [aw=avg_tot], fe
	predict viol_part if e(sample), e
	replace viol_part = viol_part * 100 //to convert to percentage points (already in fraction of mean from original paper)
	qui xtreg tmean prcp i.month juris#c.doy juris#c.doy2 juris#c.doy3 juris#c.doy4 if viol ~= . [aw=avg_tot], fe
	predict temp_part if e(sample), e
	replace temp_part = temp_part*5/9 // convert from F to C
	egen state_month_year = group(state_fips month year) // generating a state-by-year variable for use as a blocking var in the block bootstrap for the watercolor reg
	outsheet viol_part temp_part state_fips month year state_month_year using "$wd\output\DataForFigure2\JACOB.csv", comma replace

clear


***** Larrick and Timmerman 2011
	use "$lt\Larrick_Timmerman_2011.dta"

	areg hbp c_Temperature i.year if retaliation_possible == 1, absorb(ParkID) cluster(gameid)
	local eff = _b[c_Temperature]
	loneway c_Temperature ParkID if e(sample)
	local sd = r(sd_w)
	summ hbp if e(sample)
	di `eff'/r(mean)*`sd'*100

	* For Fig 2 plot
	outsheet hbp_r temp_r using "$wd\output\DataForFigure2\LARRICK.csv", comma replace

clear
	
	

***** Levy et al 2005 

	use "$levy\Levy_etal_2005.dta"
	tsset seqv year 

	replace yr_tot = . if yr_tot == -9999
	replace yr_tot = -yr_tot
	replace outbr1 = . if outbr1 == -9999
	replace outbr2 = . if outbr2 == -9999
	replace outbr3 = . if outbr3 == -9999

	gen any_outbr = (outbr1 == 1 | outbr2 == 1 | outbr3 == 1)
	xi i.year, prefix(_y_)

	//so units in % of mean and generate residualized variables:
	sum any_outbr
	loc mean_any_outbr = r(mean)
	replace any_outbr = any_outbr/`mean_any_outbr'*100 

	qui areg any_outbr _y_* if yr_tot ~= ., absorb(seqv)
	predict conflicts_xb if e(sample), xbd
	gen conflicts_r = any_outbr - conflicts_xb

	qui areg yr_tot _y_* if any_outbr ~= ., absorb(seqv)
	predict yr_tot_xb if e(sample), xbd
	gen yr_tot_r = yr_tot - yr_tot_xb

	egen country_year = group(unmj year)

	* Standardized effects
	* with clustering at the grid and country-year level
		xtivreg2 any_outbr yr_tot _y_*, fe ivar(seqv) cluster(seqv country_year) small 
		local eff = _b[yr_tot]
		loneway yr_tot seqv if e(sample)
		local sd = r(sd_w)
		summ any_outbr if e(sample)
		di `eff'/r(mean)*`sd'*100

	* Export residualized variables for figure 2:
	keep conflicts_r yr_tot_r country_year seqv
	duplicates drop
	drop if conflicts_r ==.
	drop if yr_tot_r ==.
	drop if country_year ==.
	drop if seqv ==.

	outsheet conflicts_r yr_tot_r country_year using "$wd\output\DataForFigure2\LEVY.csv", comma replace

clear


***** Buhaug Theisen Holtermann 2011 

	* Standardized effects, clustering on cell and country-year
	use "$bht\Buhaug_Theisen_Holtermann_2011.dta"
	egen cty_yr = group(gwcode year)
	xi i.year
	xtivreg2 onsetx spi6dum _I*, fe ivar(poly_id) cluster(poly_id cty_yr) small //using xtivreg for multi-way clustering. 
	local eff = _b[spi6dum]
	local se = _se[spi6dum]
	loneway spi6dum poly_id
	local sd = r(sd_w)
	summ onsetx
	di `eff'/r(mean)*`sd'*100  //standardized effect

clear


***** Ranson 2012

	use "$ra\Ranson_2012.dta"
	egen state_year = group(state year) 
	replace cmcy_maxtemp = cmcy_maxtemp*5/9 //convert to C

	//county-by-month and county-by-year effects have already been removed

	* Rape
	cgmreg cmcy_rate_rape cmcy_maxtemp cmcy_precip cmcy_lag_maxtemp cmcy_lag_precip [aweight = pop], cl(county state_year)
	local eff = _b[cmcy_maxtemp]
	keep if e(sample)
	sum cmcy_maxtemp [aweight = pop] if e(sample) //FE have already been partialed out
	local sd = r(sd)
	sum rate_rape [aweight = pop] if e(sample)
	di `eff'/r(mean)*`sd'*100

	* Murder
	cgmreg cmcy_rate_murder cmcy_maxtemp cmcy_precip cmcy_lag_maxtemp cmcy_lag_precip [aweight = pop], cl(county state_year)
	local eff = _b[cmcy_maxtemp]
	keep if e(sample)
	sum rate_murder [aweight = pop] 
	di `eff'/r(mean)*`sd'*100

	* Assault
	cgmreg cmcy_rate_assaultaggr cmcy_maxtemp cmcy_precip cmcy_lag_maxtemp cmcy_lag_precip [aweight = pop], cl(county state_year)
	local eff = _b[cmcy_maxtemp]
	keep if e(sample)
	sum rate_assaultaggr [aweight = pop] if e(sample)
	di `eff'/r(mean)*`sd'*100

	* For plots
	qui reg cmcy_rate_rape  cmcy_precip cmcy_lag_maxtemp cmcy_lag_precip [aweight = pop] 
	predict rate_rape_resid if e(sample), resid
	qui reg  cmcy_maxtemp cmcy_precip cmcy_lag_maxtemp cmcy_lag_precip [aweight = pop]
	predict temp_resid if e(sample), resid
	sum rate_rape [aweight = pop] if e(sample)
	loc mean_rate_rape = r(mean) 
	replace rate_rape_resid = rate_rape_resid*100/`mean_rate_rape' // so units are in % of mean
	outsheet rate_rape_resid temp_resid state_year using "$wd\output\DataForFigure2\RANSON.csv", comma replace

clear


***** Fjelde and von Uexkull 2012

	use "$fvu\Fjelde_vonUexkull_2012.dta"
	
	* OLS version of their Model 5 in Table 1, without controls, with year FE, and clustering at the commune level
	* rainfall data already standardized, but then censored at 0, so actual sd < 1

	xtreg xns23 negSDdetail_rev_max i.year, fe i(adm1_code) cl(adm1_code)
	loneway negSDdetail_rev_max adm1_code
	loc sd = r(sd_w)

	summ xns23
	di _b[negSDdetail_rev_max]/r(mean)*`sd'*100  //standardized effect.  

clear


log close
