clear all
set processor 2
cd "F:\01_My_work_pre\00_Project\17_GME\01_Linear\Code"
set obs 100000
gen x1 = rnormal()
gen x2 = rnormal()
gen e = rnormal()
gen y = x1 + x2 + e
program gme, plugin
tempname ze
tempname zv
tempname beta
tempname vcovm
tempname std
tempname std1
qui sum y
mat `ze' = (-3*`r(sd)'*`r(sd)'\ 0 \ 3*`r(sd)'*`r(sd)')
mat `zv' = (-5,0,5\-5,0,5)
mat `beta' = J(2,1,0)
mat list `beta'
mat `vcovm' = J(2,2,0)
mat list `vcovm'
plugin call gme y x1 x2, `zv' `ze' `beta' `vcovm'
mat list `beta'
mat list `vcovm'
di me
mat `std1' = vecdiag(`vcovm')
matewmf `std1' `std', function(sqrt)
mat list `std'
reg y x1 x2, noconstant
