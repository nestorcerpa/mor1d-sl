## Software : mor1d-sl


## Licence 

MIT Licence

Copyright (c) 2019  
    Nestor Cerpa [nestor.cerpa@gm.univ-montp2.fr]  
	Geosciences Montpellier  
	University of Montpellier  
	Place Eugene Bataillon  
	34095 Montpellier Cedex05  
	France  
    David Rees Jones [david.reesjones@earth.ox.ac.uk]   
	University of Oxford  
	Earth Sciences  
	South Parks Road  
	Oxford, OX1 3AN  
	United Kingdom  
    Richard F. Katz [richard.katz@earth.ox.ac.uk]  
	University of Oxford  
	Earth Sciences  
	South Parks Road  
	Oxford, OX1 3AN  
	United Kingdom  

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## REFERENCE: 
Cerpa, N. G., Rees Jones, D. W., and Katz, R. F.:
Consequences of glacial cycles for magmatism and carbon transport at mid-ocean ridges 
(Please cite final revised manuscript when available)

## CODE SUMMARY: 
* All functions needed to compute the 1-d melting column models are provided in src/
* Functions mean_analytical_xxx.m solve eqs. (A.1) and flucuations_xxx.m solve eqs. (B.1)   
* The model parameters are defined in input_parameters.m
* Data files used for plotting Figures 5, 8 and 9 are provided 

## INSTRUCTIONS:
Follow the instructions below to generate the figures in the EPSL paper:  
1. To generate Figures 2 and 3: Run mor1d.m  
2. To generate Figures 4, 6(a,b,c) and 7(a,b,c): Run mor1d_wetdry_fluxes.m  
3. To generate Figures 5(a,b), 6(d,e) and 7(d,e): Run mor1d_admittance_plot.m (to re-calculate data used by the script : run mor1d_admittance_compute.m)   
4. To generate Figures 5(c,d): Run mor1d_admittance_withQ_plot.m (to re-calculate data used by the script: run mor1d_admittance_withQ_compute.m)  
5. To generate Figures 8 and 9 from the EPSL paper : Run mor1d_slrecord_plot.m (to re-calculate data used by the scipt : run mor1d_slrecord_compute.m)  
  
n.b.: To generate Figures S1 to S5 in the Supplementary Material change parameters in input_parameters.m accordingly and run mor1d.m; To generate Figure S6 run mor1d_asymptotics_plot.m  

## Additional notes 

* This distribution of mor1d-sl was developed and tested on MATLAB2018b. 
Running these routines on different versions of MATLAB may lead to compatibility issues.
