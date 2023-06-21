# RSPA-2023-0195
<b>A Semi-Analytic Method for the Computation of the Effective Properties</b>

This is a repository of the codes used for this manuscript.

<b>Example 1: Disk of area fraction $f< f_crit$ in square/hexagonal array</b>

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
