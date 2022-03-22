
****************************************************
*** Linear-mixed models in economics evaluations ***
****************************************************


/*
Notes:
	Stata code supplement to 2022 Health Economics paper by Gabrio et al. 
	"Linear mixed models to handle missing at random data in trial-based economic evaluations" (Gabrio et al., Health Economics 2022)
	Author: Baptiste Leurent
	Date: March 2022
	Stata version: 16
	Dataset used: SADD trial (ISRCTN88882979)
*/


** Open data
	use "data_long_SADD.dta", clear
	gen trtp = trt*(time>1) //Treatment arm indicator, for follow-up visits only
	
** LMM for utilities
	mixed u i.time i.time#i.trtp || id:, res(unstructured, t(time))
		// || id: = random effect for each participant
		// res(unstructured, t(time)) = unstructured covariance, for each timepoint
		// i.time#i.trtp = interactions between time and treatment
	
	*Incremental utility at each time point:
		lincom 2.time#1.trtp  //Treatment effect (trtp=1) at time 2 
		lincom 3.time#1.trtp  //Treatment effect at time 3
	*Incremental QALYs
		lincom 0.375*2.time#1.trtp + 0.25*3.time#1.trtp
	
	*Utility per arm at each time-point
		margins i.time#i.trtp
	*QALYs per arm
		margins i.time#i.trtp, post
		lincom 0.125*1.time#0.trtp + 0.375*2.time#0.trtp + 0.25*3.time#0.trtp
		lincom 0.125*1.time#0.trtp + 0.375*2.time#1.trtp + 0.25*3.time#1.trtp
	
** LMM for costs
		mixed c i.time i.time#trtp || id:, res(unstructured, t(time))
		
		*Incremental cost at each time point:
			lincom 2.time#1.trtp  //Treatment effect at time 2
			lincom 3.time#1.trtp  //Treatment effect at time 3
		*Incremental total cost
			lincom 2.time#1.trt + 3.time#1.trt  
	
		*Cost per arm at each time-point
			margins i.time#trtp
		*Total cost per arm
			margins i.time#trtp, post
			lincom 2.time#0.trt + 3.time#0.trt
			lincom 2.time#1.trt + 3.time#1.trt
	