# RSPA-2023-0195
<b>A Semi-Analytic Method for the Computation of the Effective Properties of Composites of Two Isotropic Constituents</b>

This is a repository of the codes used for this manuscript.

<b>Example 1: Disk of area fraction f< f_c in square/hexagonal array</b>

code: spectral_disk.f90

<b>Usage:</b> 

./spectral_disk [geom] [f]  [Ng] [Ndpt]  [Nmom]
    
  geom: S/T (square/hexagonal) 
  
  f: area fraction
  
  Ng : Number of Green function coefficients
  
  Ndpt : Number of boundary points per unit length
  
  Nmom: Number of moments (even integer)
  
<b>Example:</b> ./spectral_disk S 0.7 6 500 10

</b>Output: the moments mu[k] k=1..Nmom</b>

 Npt =         1482  boundary points <br>
 series coefficients for sigma*<br>
 
 mu[1]:=  0.70000000000001183      ;<br>
 mu[2]:= -0.10499999999999998      ;<br>
 mu[3]:=   2.8735366212673228E-002 ;<br>
 mu[4]:=  -1.2750792970138583E-002 ;<br>
 mu[5]:=   6.9749243309943998E-003 ;<br>
 mu[6]:=  -4.1692813060702294E-003 ;<br>
 mu[7]:=   2.6174142590385373E-003 ;<br>
 mu[8]:=  -1.6957384466021780E-003 ;<br>
 mu[9]:=   1.1221268517625290E-003 ;<br>
 mu[10]:=  -7.5341424121961893E-004 ;<br>
  
 Exec time  =   0.26382304599974304      seconds<br>

 To build the Pade approximants in the variable s=(epsilon1-epsilon2)/epsilon2, run:
 
 maple pade_spect_disk.mw
 
</b>Output: PA values at s=s0, s=-1 and s-> infinity </b>

"sub/diagonal PA at s=", 50 <br>
                                         2, -0.02941176470588184, 5.117647058823661<br>
                                          4, -1.559076082142545, 6.211205816742471<br>
                                          6, 11.55626291207727, 6.348268279825770<br>
                                          8, 6.486991598507176, 6.353978155550265<br>
                                          10, 6.356935150708070, 6.35409689070327<br>
                                          
"sub/diagonal PA at s=-1"<br>
                                         2, 0.5882352941176430, 0.1764705882352827<br>
                                         4, 0.1576774758102364, 0.1388843633851914<br>
                                         6, 0.1368071384146264, 0.1347299134440655<br>
                                         8, 0.1346372886161442, 0.1345446637882222<br>
                                         10, 0.1345425709654695, 0.1345404781427172<br>
 "diagonal PA at s=infinity"<br>

                                                    2, 5.666666666666825<br>
                                                    4, 7.200234609755841<br>
                                                    6, 7.422256679584068<br>
                                                    8, 7.432476114950040<br>
                                                   10, 7.432707344320310<br>
