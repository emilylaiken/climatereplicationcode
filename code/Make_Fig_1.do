
* This file creates Figure 1 in Hsiang Burke Miguel 2013

clear all

* replace the following line with the directory where the replication data was unzipped

global wd "/Users/solhsiang/Dropbox/SHARED_FOLDERS/Marshall-Sol/drafts/Science_review_v4/replication/"

cd $wd
insheet using data/Figure1_data.csv

encode continent, gen (continent_ID)
sort continent_ID start_year
gen N = _n

sum N if continent == "AF"
loc N_AF = r(mean)
sum N if continent == "AM"
loc N_AM = r(mean)
sum N if continent == "EA"
loc N_EA = r(mean)
sum N if continent == "GB"
loc N_GB = r(mean)

gen log10_wavelength = log10(wavelength_in_years)
gen mid_year = (end_year - start_year)/2+start_year
gen log10_start_year_ybp = -log10(2012 - start_year )
gen log10_end_year_ybp = -log10(2012- end_year)
gen log10_mid_year_ybp = -log10(2012- mid_year)
gen mid_log_years = (log10_end_year_ybp-log10_start_year_ybp)/2+log10_start_year_ybp  


loc ybc8000 = -log10(2012+8000)
loc ybc4000 = -log10(2012+4000)

loc year_list "0 1000 1500 1800 1900 1950 1980 2000 2010"


foreach i of loc year_list{
	loc y`i'= -log10(2012-`i')
}



tw (rbar log10_start_year_ybp log10_end_year_ybp N if continent == "AF",  horizontal color(gs12))(rbar log10_start_year_ybp log10_end_year_ybp N if continent == "AM", horizontal color(gs7))(rbar log10_start_year_ybp log10_end_year_ybp N if continent == "EA", horizontal color(gs3))(rbar log10_start_year_ybp log10_end_year_ybp N if continent == "GB", horizontal color(black)),  xlab(`ybc8000' "8000 BCE" `y0' "0" `y1000' "1000"  `y1800' "1800" `y1950' "1950" `y2000' "2000" `y2010' "2010") xline(`ybc8000'  `y0'  `y1000'  `y1500' `y1800' `y1950'    `y2000'  `y2010' , lcolor(gs11)) xtit("Years in study (log scale in YBP)") ytit("Region") legend(off)  ylab(`N_AF' "Africa" `N_AM' "Americas" `N_EA' "Eurasia" `N_GB' "Global", angle(0))
graph export "output/Figure1a.pdf", as(pdf) replace

gen spatial_wavelength = .
replace spatial_wavelength = -1 if spatial_scale == "site"
replace spatial_wavelength = 0 if spatial_scale == "municipal"
replace spatial_wavelength = 1 if spatial_scale == "pixel"
replace spatial_wavelength = 2 if spatial_scale == "state"
replace spatial_wavelength = 3 if spatial_scale == "country"
replace spatial_wavelength = 4 if spatial_scale == "region"
replace spatial_wavelength = 5 if spatial_scale == "global"
replace spatial_wavelength = 4.5 if study == "Zhang et al 2011"

egen time_space = group( spatial_wavelength log10_wavelength )
bysort time_space : egen N_time_space = count(time_space )
sort start_year

tw sc spatial_wavelength log10_wavelength [aweight = N_time_space], xlabel(-3.9425051 "hour" -2.5622929 "day" -1.69897 "week" -1.079 "month" 0 "year" 1 "decade" 2 "century" 3 "millennia" , angle(90) tlcolor(gs8))  ylabel(-1 "site" 0 "municipal" 1 "pixel" 2 "state" 3 "country" 4 "region" 5 "global" , angle(0) tlcolor(gs8)) m(o) msize(medsmall) mlwidth(vthin) mfcolor(emidblue ) mlcolor(black) xtit(" " "Wavelength of climate signal (log scale)") ytit("Spatial scale of dependent variable" " ") yscale(lstyle(none)) xscale(r(3.1) lstyle(none))  xsize(1) ysize(1) plotregion(style(none))
graph export "output/Figure1b.pdf", replace


