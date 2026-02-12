
Program                         InterdigitatedFilter ;

  {Thie program creates round rod interdigitated microwave filters with 
  equal diameter rods.  Input and output are rod coupled. 
  
  Interdigitated filter from "Design of Wide-Band (and Narrow-Band)
  Band-Pass Microwave Filters on the Insertion Loss Basis,"  Matthaei
  IRE Transactions on Microwave Theory and Techniques, Nov 1960 pp 580-593}
  
  {Created 1997 but not released}
  {Revised 10/24/2018, not released}
  {Verion 1-0: 02/01/2021, Added loss and rod length estimators, data for
     QUCS simulation.  first released version}
  
  { Copyright (c) 2021  Dennis G. Sweeney  wa4lpr@arrl.net

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See <https://www.gnu.org/licenses/ for details of the 
    GNU General Public License.
}

{Complied with the Free Pascal Complier (FPC), see www.freepascal.org.  
 Versions of FPC are available for just about any processor or operating system
 you can think of.  The Ingrated Development Environment (IDE) used is Geany,
 see www.geany.org.  Geany works with C, Python, Pascal and a number 
 of other programming languages.  It is availible for Linux, Mac OSX, and 
 Winders.  Both FPC and Geany are available for no charge under GPL.

All calculations done in cm , Hz.
}
 
  
uses crt, math ;

type
      vector = array[0..20] of real ;

var                    	N				: integer ;	{number of filter elements}
						Ans		        : integer ;	{Answer}						
						j				: integer ; {counter index}	
						
						b				: real ;	{ground plane spacing}
						dia				: real ;	{rod diameter}
						bprime			: real ;	{outer dia of equivalent coax transformed from slab line}
						rprime			: real ;	{normalized inner radius of equivalent co-ax}
						lamda0			: real ;	{center frequecy wavelength}
						e				: real ; 	{Filter fipple} 
						Q, P        	: real ;	{Intermediate values, paper Table I}
						F		        : real ;	{Filter center frequency}
						BW          	: real ;	{Filter bandwidth}
						Qf              : real ;	{Filter Q = F/BW}
						
                        g				: vector ;	{protype g values} 
                        Yoo, Yoe   		: vector ;	{Even and Odd mode Admitances}
						Zoo, Zoe		: vector ;	{Even and Odd mode Impedances}
						m, K, NN        : vector ;  {coupling coefficient, impedance inverters & normalized coupling coefficient,  }
													{intermediate value, paper Table 1}
						Zeven, Zodd     : real ;	{Zoe/Zoo for fixed d/h}
						Zoeprime, Zooprime : real ;	{Ze/Zo for modeling}
						Zs              : real; 	{modeling Zo}
						Loss			: real ;    {estimated filter loss}
						Zo              : real ;	{Input/Output Impedance}
						ZoSlab			: real ;	{Zo of slab line}
						ZoCoax			: real ;	{Zo of z = tan(w) transformer coax}
						dh				: real ;	{filter rod dia / ground plan spacing}
						sh				: real ;	{s/h resonator spacing }
						cc				: real ;	{rod center to center distance}
						l				: real ;	{resontator length}
						resonatorQ		: real ;	{resonator Q}
						gap				: real ;	{distance from rod end to filter wall}
						spacing			: real ;	{chosen gap spacing}
						Cp, Cf, Ct		: real ;	{rod-wall, fringing, total Cap}
						SuMg			: real ;	{sum of prototype g values for loss}
						Theta1, Ya		: real ; 
						s              	: real ;  	{Impedance scale factor}
						

{Chebyshev       returns  g: array of low pass prototype values for
                          Chebyshev response, g is a global variable because
                          you can't pass an array as a varible in a proceedure.
                 input    N: number of elements
                          ripple: passband ripple in dB}

Procedure Chebyshev(N:integer ; ripple : real) ;

var     A, B       : array [0..20] of real ;
        C          : real ;
        Beta, Gamma: real ;
        k          : integer ;

begin
  C := 2.0 * ripple / 17.37 ;
  Beta  := ln((exp(C) + 1.0) / (exp(C) - 1.0)) ;
  Gamma := 0.5 * (exp(Beta / (2.0 * N)) - exp(-Beta / (2.0 * N))) ;
  for k := 1 to N do
  begin
    A[k] := sin(0.5 * (2.0 * k - 1.0) * pi / N) ;
    B[k] := sqr(Gamma) + sqr(sin (k * pi / N)) ;
  end ;
  g[0] := 1.0 ;
  g[1] := 2.0 * A[1] / Gamma ;
  for k := 2 to N do
       g[k] := 4.0 * A[k-1] * A[k] / (B[k-1] * G[k-1]) ;
  If Odd(N) then
    g[N + 1] := 1.0
  else
    g[N + 1] := sqr((exp(Beta / 2.0) + 1.0) / (exp(Beta / 2.0) - 1.0)) ;
end ;

{TheedBBW calculates the ratio of 3dB BW to ripple BW for a Chebyshev filter
	Input : 	ripple: passband ripple in dB
				N: 		number of elements
    Returns :ratio of 3dB BW / Ripple BW }

function ThreedBBW(N : integer ; ripple : real) : real ;

var		delta : real ;

begin
	delta := sqrt(power(10.0,(0.1*ripple))-1.0) ;
	ThreedBBW := cosh((1.0/N)*arcosh(1.0/delta)) ;
end ;	

{RippleRL returns the in band return loss for a given ripple for a 
Chebyshev filter
	Input:	ripple: passband ripple in dB
	Output: Return Loss in dB}	

function RippleRL(ripple : real) : real ;
begin
	RippleRL := -10.0 * log10(power(10.0,ripple/10.0) - 1.0);	
end ;		

{Butterworth     returns  g: array of low pass prototype values for
                          Butterworth response, g is a global variable because
                          you can't pass an array as a varible in a proceedure.
                 input    N: number of elements }


procedure Butterworth (N : integer) ;

var       r            :integer ;
          Theta        : real ;

Begin
   g[0] := 1.0 ;
   g[N + 1] := 1.0 ;
   For r := 1 to N Do
   Begin
      Theta := (2 * r - 1) * Pi / (2 * N ) ;
      g[r] := 2 * Sin(Theta) ;
   End ;
End ;


{Slab line impedance (in air) 
	Input: 	spacing: slab line ground plane spacing
			dia: rod diameter
	Outputs: characteristic Zo

Rosloniec, Stanislaw, "Algorithms for 
Computer-Aided Design of Linear Microwave Circuits," Artech House
Boston, 1990, section 6.3}

function slabZo(spacing,dia : real) : real ;
var

	R, X, Y : real ;
	
begin	
	R := 0.25 *  pi * dia/spacing ;
	X := 1 + 2.0 * sqr(sinh(R)) ;
	Y := 1.0 - 2.0 * sqr(sin(R)) ;
	slabZo := 59.952 * ln((sqrt(X) + sqrt(Y))/sqrt(X-Y)) - sqr(R)*sqr(R)/30.0 + 0.014 * sqr(sqr(R)*sqr(R)) ;
end ;

{Zo of coax line (in air)
	Input: 	outerdia: outer diameter of coax line
			innnerdia: diameter of center conductor
	Outputs: characteristic Zo}

function coaxZo(outerdia,innerdia :real) : real ;
begin
	coaxZo := 60.0 * ln(outerdia/innerdia) ;
end ;

{endCap calculates rod capacitive loading.  This is for a coax line. The 
slab line is transformed with z = tan(w) into a coax approximation.
	Cf:	fringing capacitance	Cp: parallel plate capacitance from rod to wall
 	Input: 	outerdia: coax outer diameter
 			innerdia : coax inner diameter
			gap : distance from the end of the rod to the wall
 	Ouput: Cf + Cp in pf (mult by 1e-12 for F) }
 	
{Nicholson, B.F., "The Resonant Frequency of Interdigital 
Filter Elements," IEEE Transaction on Microwave Theory and Techniques,
Vol. MTT-14, No. 5, May 1966, 250-251.}

{A 5th order polynomial curve fit is used to approximate Nicholson's 
fringing capacitance.  The data was obtained by measuring the plot in the 
above paper with digital calipers and then doing a polynomial curve fit.  
The curve fit coefficients were obtained by using Octave's polyfit}
	
function FringCap(outerdia,innerdia : real) : real ;
var 
		w :	 		real ;
begin
	w := innerdia/outerdia ;
	
	FringCap := pi * outerdia * (0.1484*w*sqr(sqr(w)) - 0.13989*sqr(sqr(w)) + 0.12623*w*sqr(w) - 0.05089*sqr(w) + 0.089164*w - 2.8858e-5) ;
end ;

{EndCap calculates rod parallel capacitance between rod end and filter wall.
This is for a coax line. The slab line is transformed with z = tan(w) into
a coax approximation.
 	Input: 	innerdia : coax inner diameter
			gap: distance between rod end and filter wall
			Since inner/gap ratio is used, must be the same units.			
 	Ouput: EncCap in pf (mult by 1e-12 for F) }
 	
{Nicholson, B.F., "The Resonant Frequency of Interdigital 
Filter Elements," IEEE Transaction on Microwave Theory and Techniques,
Vol. MTT-14, No. 5, May 1966, 250-251.}

function EndCap(innerdia, gap : real) : real ;
begin	
	EndCap := 0.0885*pi*sqr(innerdia)/(4*gap) ;  {rod parallel plate cap}
end ;

{Length	: 1/4 wave resonator Length loaded with an end capacity C
	Input: 	Zo:	resonator impedance
			fo:	operating frequency
 		 	C:	end capactity in pF
	Returns : L: resonator length in cm }
	 
{Nicholson, B.F., "The Resonant Frequency of Interdigital 
Filter Elements," IEEE Transaction on Microwave Theory and Techniques,
Vol. MTT-14, No. 5, May 1966, 250-251.}

function Length(Zo, fo, C : real ) : real ;

begin 
	Length := (100 * 299.792e6 / (2 *pi * fo)) * arctan(1/(2*pi*fo*Zo*(C)*1e-12)) ;
end ;
					
{rodQ calculates the Q of a 1/4 wavelengh rod }

{Nicholson, B.F., "Dissipation Loss in Interdigital and Combline Filters,"
Electronics Letters, March 1966, Vol. 2, No. 3, pp 90-91.}
{
Input:	r: rprime: normalized inner rod,innerdia/outerdia for equvalent coax
		b: ground plane spacing cm
		l: finger length cm
		fo: frequecy in MHz
Ouput:	resonator Q

Note: The "FudgeFactor" is suggested by Nicholson to accomidate real world Q.
Nicholson suggests 0.8, but 0.6 seems more consistent with measurement on 
CU/Brass and AL/Brass filters.  Your milage may (will!) vary, but assume
that Q is always less than theory.		
}

function rodQ(r, b, l, lamda0, fo :real) : real ;
var
	mu, ro :		real ;
	FudgeFactor :	real ;
begin

	mu := 1.257e-8 ;  		{h/cm}
	{ro := 1.724e-6 ;}		{ohm-cm for copper}
	{ro := 1.591e-6 ;}		{ohm-cm for silver}
	{ro := 4.066-6 ;}		{ohm-cm for 6060-T6 AL}
	ro := 6.897e-6 ;		{ohm-cm for yellow brass}
	{http://eddy-current.com/conductivity-of-metals-sorted-by-resistivity/}
	{Pick your favorite flavor!}
	
	FudgeFactor := 0.80 ;
	
	rodQ := (FudgeFactor) * 240.0 * sqr(pi) * ln(1/r)/((0.5*pi*lamda0/b)*
			sqrt(pi*fo*mu*ro)*(1+1/r) + 
			8.0*sqr(0.25*lamda0/l)*sqrt(pi*fo*mu*ro)*ln(1/r)) ;
end ;	

{Couple computes s/h from d/h and coupling coefficient K for
 1/4 wavelength round rods

			returns:	y = s/h:  rod edge to edge spacing
						Zeven, Zodd impedance of coupled rods for
						given d/h and coupling coefficient.
			input:		x = d/h	
						K = coupling coefficient 
}				
{Couple uses Even and Odd mode impedances for coupled slabline.

 Stanislaw Rosloniec, "Alogrithms for Computer-Aided Design of Linear 
 Microwave Circuits," Artech House, 1990, pp. 203-208.
 
 Stanislaw Rosloniec, "An Improved Algorithm for the Computer-Aided 
 Design of Coupled Slab Lines," IEEE Transaction on MTT, Vol 37, No. 1, 
 January 1989, pp 258-261. 
}

procedure Couple(var x, y, K, Zeven, Zodd : real) ;


var   y1, y2, yn     : real ;
      ff1, ff2, ffn   : real ;
      n              : integer ;

  function f1(x : real) : real ;
  begin
    f1 := x * (1.0 + exp(16.0 * x - 18.272)) / sqrt(5.905 - sqr(x) * sqr(x)) ;
  end ;

  function f2(x,y : real) : real ;
  var    y2, y3   : real ;
  begin
    y2 := sqr(y) ; 
    y3 := y2 * y ;

    if y < 0.9 then
      f2 := (-0.8107 * y3 + 1.3401 * y2 - 0.6929 * y + 1.0892 +
                                              0.014002 / y - 0.000636 / y2) -
           x * (0.11 - 0.83 * y + 1.64 * y2 - y3) +

          -0.3345 * exp(-13.0 * x) * exp(- 7.01 * y + 10.24 * y2 - 27.58 * y3)
    else
      f2 := 1.0 + 0.004 * exp (0.9 - y) ;
  end ;

  function f3(x, y : real) : real ;
  
  begin
     
     f3 := tanh(0.5*pi*(x+y)) ;
    
  end ;

  function f4(x, y : real) : real ;
  var   y2 , y3    : real ;
  begin
    y2 := sqr(y) ;  y3 := y2 * y ;
    if y < 0.9 then
       f4 := (1.0 + 0.01 * (-0.0726 - 0.2145 / y + 0.222573 / y2 - 0.012823 / y3))
          - x * 0.01 * (-0.26 + 0.6866 / y + 0.0831 / y2 - 0.0076 / y3) +
          (-0.1098 + 1.2138 * x  - 2.2535 * sqr(x) + 1.1313 * x * sqr(x)) *
          (-0.019 - 0.016 / y + 0.0362 / y2 - 0.00243 / y3)
    else
      f4 := 1.0 ;
  end ;

  function Zoe(x,y : real) : real ;
  begin
    Zoe := 59.952 * ln(0.523962 / (f1(x) * f2(x,y) * f3(x,y))) ;
  end ;

  function Zoo(x,y : real) : real ;
  begin
    Zoo := 59.952 * ln(0.523962 * f3(x,y) / (f1(x) * f4(x,y))) ;
  end ;
  
  {Secant method back solve will drive this function to 0, resulting 
  in a value of y (s/h) that gives the specified value of coupling}

  function f(x,y,k : real) : real ; 
  begin
    f := (Zoe(x,y) - Zoo(x,y)) / (Zoe(x,y) + Zoo(x,y)) - k ;
  end ;

{Secant method to back solve. This is a quasi Newtons Method algorithm.
* You have x=d/h which is fixed and k. You have k = f(d/h,s/h) and you need
* y=s/h for a given k value.  This must be solved iteratively.   
* CAUTION: The secant method is supposted to converge faster 
* than simple bracketing methods but it may not converge for some values.  
* If convergance fails, you will get a "failed to converge" error message.  
* The resulting values may still be usable but with degraded accuracy.
* 
* 
* See Johnson, Lee W. and R. Dean Riess, "Numerical Analysis," 2nd Ed, 
* Addison-Wesley Publishing, Massachusetts, 1982,  p166-168 
* or 
* Flannery, Teukolsky, Vetterling, "Numerical Recipes in Pascal," 
* Cambridge Univeristy Press, Cabriadge, 1989, University Press, 
* Cambridge, 1989, Section 9.2. }

begin
	{x is d/h and fixed, y is s/h}
	
	y1 := 0.1 ;		{low root estimate}
	y2 := 3 ;		{high root estimate }
	n := 0 ; 		{iteration counter}

	repeat
		n := n + 1 ;
		ff1 := f(x,y1,K) ;  {f(y1) }
		ff2 := f(x,y2,K) ;	{f(y2) }
		yn := (y1 * ff2 - y2 * ff1) / (ff2 - ff1) ; {y3 = (y1*f(y2)-y2*f(y1))/(f(y2)-f(y1))}
		ffn := f(x,yn,K) ;   {f(y3) }

		if (ff1 * ffn) < 0 then {f(y3) * f(y1) < 0}
		begin
			y2 := yn ;
			ff2 := ffn ;
	end
	else
	begin
		y1 := yn ;
		ff1 := ffn ;
	end ;
until (abs(ff1*ffn) < 1e-9) or (n > 600) ;

if n > 600 then 
	begin
			writeln('Error: failed to converge') ;
			writeln('n = ',n,' ff1*ffn = ',abs(ff1*ffn)) ;
	end ;
			
	y := yn ;

	Zeven := Zoe(x,y) ;  Zodd := Zoo(x,y) ;
end ;

{fx function 
* x: gap between rod end and filter wall
* calculate Cp for value of x
* Spacing, dia, bprime, ZoCoax, F are global variables}

function fx(x:real) : real ;
var Cp, Ct, l : real ;
begin
	Cp := EndCap(dia,x) ;
    Ct := FringCap(bprime, dia) + Cp ;
    l := Length(ZoCoax, F, Ct)	;
    fx := spacing-(x + l) ;
end ; 
   
{Bisection back solve, bisection is a bracking algorithm.  It is slower
than the Secant method above but it should always converge.  It will return
a failure to converge message if it doesn't work.

	Input: 	x1: solution low limit		solution must be between x1 and x2.
			x2: solution upper limit
			xacc : error tolerance
	Calls:	function fx
	Returns: value of x that forces f(x) = 0
	* 
Flannery, Teukolsky, Vetterling, "Numerical Recipes in Pascal," 
Cambridge Univeristy Press, Cabriadge, 1989, p 277. }

function rtbis(x1,x2,xacc: real): real;

var
      dx,f,fmid,xmid,rtb: real;  j: integer;
begin
      fmid := fx(x2);  
      f := fx(x1);
      if (f < 0.0) then 			{decide which end to start from}
      begin
         rtb := x1; 
         dx := x2-x1  
      end
      else 
      begin
        rtb := x2;  
        dx := x1-x2  
      end;
      
      j := 0 ;
      repeat
        dx := dx*0.5; 				{cut the interval in half}
        xmid := rtb+dx;  
        fmid := fx(xmid);
        if (fmid <= 0.0) then rtb := xmid;
        j := j + 1 ;					{count how many times through loop}
      until ((abs(dx) < xacc) or (fmid = 0.0) or (j > 40)) ;
        {if the interval is less than error tolerance or more the 40 loops, quit}
        
      if (j > 40) then 
      begin
		write('Failed to converge after 40 iterations') ;
		readln ;
	  end ;
      
	  rtbis := rtb ;
end;

BEGIN  {Main Program}

	Clrscr ;

	writeln('INTRFIL: Interdigitated Filter Program  Version 1.0') ;
	writeln('Copyright (c) 2021 Dennis G. Sweeney  WA4LPR') ;
	writeln ;
	write('Number of elements? ') ; 		readln(N) ;	
	write('Center frequency (MHz) ') ; 		readln(F) ;
	F := F * 1e6 ;
	write('Butterworth (1)  Chebychev (2) ? ') ; readln(Ans) ;

	if ans = 1 then
	begin
		write('3 dB Bandwidth (MHz) ') ; 	readln(BW) ;
		Butterworth(N) ;
	end
	else
	begin
		write('Ripple Bandwidth (MHz) ') ; 	readln(BW) ;		
		write('Bandpass ripple (DB) ') ; 	readln(e) ;
		Chebyshev(N, e) ;
		writeln('3 dB Bandwidth (MHz) = ',(BW*ThreedBBW(N,e)):4:2,'    Ripple Return Loss = ',RippleRL(e):3:1,' dB') ; 		
	end ;
	BW := BW * 1e6 ;
	
	{The Matthaei paper is for edge coupled filters.  An interdigited 
	filter is an edge coupled filter folded in half.  The folding reduces 
	the shunt element impedances in half.  The fix is to reduce the BW by 1/2}
	
	
	BW := 0.5 * BW ;							{folding the filter doubles the BW}
	Qf := F /BW ;								{Filter Q }	
							
	lamda0 := 100.0 * 299.792e6 / F ;  			{center frequency wavelength cm}
	
	write('Load impedance ') ; 				readln(Zo) ;
	Ya := 1 / Zo ;								{Input line admittance}
	write('ground plane spacing in = ') ; 	readln(b) ;
	b := b * 2.54 ;								{convert to cm}
	write('rod dia in = ') ; 				readln(dia) ;	
	dia := dia * 2.54 ;							{convert to cm}
	
	dh := dia/b ;								{d/h}
		
	bprime := 4.0 * b / pi ;					{outer diameter of equivalent co-ax }

	rprime := dia/ bprime ;						{normalized inner radius of equivalent co-ax}
	
	ZoSlab := SlabZo(b,dia)  ; 					{Zo of slab line }
	ZoCoax := coaxZo(bprime, dia) ;				{Zo of equvalent co-ax from z = tan(w) transformation}
	
	writeln ;
	
	writeln('d/h = ',dh:2:3,'   Slab Zo = ',ZoSlab:3:2,'    Co-ax Zo = ',ZoCoax:3:2) ;	


{See Table I of Design of Wide-Band (and Narrow-Band)
  Band-Pass Microwave Filters on the Insertion Loss Basis}

	K[0] := 1 / sqrt(g[0] * g[1]) ;				{Filter Impedance Inverters}
	K[n] := K[0] ;

	Theta1 := 0.5 * pi * (1 - 0.5 / Qf) ;

	Q := cot(Theta1) ;
	P := sqrt((Q * (Q * Q + 1)) / (Q + 0.5 / sqr(K[0]))) ;

	Zoe[0] := Zo * (1 + P * sin(Theta1)) ;		{Input and output coupling bars} 	
	Zoo[0] := Zo * (1 - P * sin(Theta1)) ;
	Zoo[n] := Zoo[0] ;  Zoe[n] := Zoe[0] ;		{Filter is symmetric}

	Yoo[0] := sqr(Ya) * Zoe[0] ;
	Yoe[0] := sqr(Ya) * Zoo[0] ;
	Yoo[n] := Yoo[0] ;  Yoe[n] := Yoe[0] ;

	s := Zo * sqr(P * sin(Theta1) / K[0]) ;		{Impedance scale factor}
	

	for j := 1 to N-1 do						{Calculate internal bars}
	begin
		K[j] :=  1 / sqrt(g[j] * g[j + 1]) ;	{Calculate nomalized coupling factor}
		NN[j] := sqrt(sqr(K[j]) + 0.25 * sqr(tan(Theta1))) ;
		Zoe[j] := s * (NN[j] + K[j]) ;			{Calculate Zoe/Zoo for internal bars}
		Zoo[j] := s * (NN[j] - K[j]) ;
		Yoo[j] := sqr(Ya) * Zoe[j] ;			{Calculate Yoo/Yoe for internal bars}
		Yoe[j] := sqr(Ya) * Zoo[j] ;
	end ;
	
	writeln ;
	writeln('Calculated odd/even impedances and coupling coefficients') ;
	writeln ;
	
	for j := 0 to n do							{Calculate resonator coupling coefficients}
	begin										{from Zoo and Zoe}
	
		m[j] := (Zoe[j] - Zoo[j]) / (Zoe[j] + Zoo[j]) ;
		writeln('Zoo[',j,'] = ',Zoo[j]:3:2,'    Zoe[',j,'] = ',Zoe[j]:3:2,
		'    K = ',m[j]:1:5) ;
	end ;
	
	writeln ; 
	
	
{Calculate a shunt 1/4 wave stub with 1/4 wave connecting lines filter. 
This is optional. This is the  equavalent circuit for a 
1/4 wave edge coupled filter.   If you make an interdigitated filter,
you effectively reduce the shunt stub impedances by 2.  This results 
in a filter twice as wide  as you specify.   You don't need this and 
the impdeances will be unrealizable but you can put the model in QUCS to
check your work.  See Figure 4b in the Matthaei paper for 
the equvalent stub circuit}

(*

	writeln('Z[0] = ',(1/Yoe[0]):2:3) ;					{Impedance of first shunt stub}
	for j := 1 to n do
		writeln('Z[',j,'] = ',(1 / (Yoe[j-1] + Yoe[j])):2:3) ; 
		writeln('Z[',(j+1),'] = ',(1/Yoe[n]):2:3) ;		 {Impedance of shunt stubs}

	for j := 1 to n + 1 do
		writeln('Z[',(j-1),j,'] = ',(2 / (Yoo[j-1] - Yoe[j-1])):2:3) ;
														{Imedance of connection lines}
													
*)
														
{Now calculate bar spacing. dh = d/h, d/h is fixed, m[j] is the coupling factor, returns s=s/h, this is
* the edge to edge spacing of the bars, but a new Zoo and Zoe that result from enforcing a fixed d/h.
* Use the resulting Zoo and Zoe to create an array of coupled bars in QUCS to test your results.
* All bars will be the same diameter.}

writeln ;
writeln('Odd/even impedances for fixed d/h = ',dh:2:3,' with s/h rod spacings') ;
writeln ;

	for j := 0 to n do
	begin
		Couple(dh,sh,m[j],Zeven,Zodd) ;
		cc := sh * b + dia ;				{center to center: s/h * grnd plane spacing + rod dia}
		writeln('Zoo = ',Zodd:3:2,'   Zoe = ',Zeven:3:2,'   s/h[',j,j+1,'] = ',sh:1:4,'    c - c = ',cc/2.54:3:3,' in') ;
	end ;
	
writeln ;

{Estimate resonator length}


Cf := FringCap(bprime, dia) ;	{bar end fringing capacitance}
l := Length(ZoCoax,F,Cf) ;		{rod length with only fringing cap}
gap := 0.25*lamda0 - l ;		{estimated gap if filter width is lamda0/4} 
Cp := EndCap(dia,gap) ;			{parallel plate C due to end gap}
Ct := Cf + Cp ;					{total capacitive loading}
l := Length(ZoCoax,F,Ct) ;		{rod length with total capacitance}

writeln('Parameters for filter cavity length = wavelength/4') ;
writeln ;
writeln('Ct = ',Ct:2:4,' pf   Cf = ',Cf:2:4,' pf   Cp = ',Cp:2:4,' pf') ;
writeln('l = ',l/2.54:2:3,' in    gap = ',gap/2.54:2:3,' in    lamda/4 = ',0.25*lamda0/2.54:2:3,' in') ;	
writeln ;

{Spacing can be adjusted approximitely +/- 10% from lamda0/4}

writeln('Cavity spacing width can be adjusted +/- 10% of wavelength/4') ;
write('Spacing = (in) ') ; readln(spacing) ;
spacing := spacing * 2.54 ;	

If (spacing <= 0.225*lamda0) then
	begin
		writeln('Spacing  <  0.9 X lamda/4, too small') ;
		write('New spacing = ') ; readln(spacing) ;	
	 	spacing := spacing * 2.54 ;		
	end
else
If (spacing >= 0.275*lamda0) then
	begin
		writeln('Spacing > 1.1 X lamda/4, too large') ;
		write('New spacing = ') ; readln(spacing) ;
		spacing := spacing * 2.54 ;
	end ;

gap := rtbis(0.001,0.35,0.0002) ; 	{Find a new gap beteen the rod end
									and the filter cavity wall.}
writeln ;

Cf := FringCap(bprime, dia) ;
Cp := EndCap(dia,gap) ;
Ct := Cf + Cp ;
l := Length(ZoCoax,F,Ct)	;

writeln('Parameters for Filter cavity ') ;
writeln ;
writeln('Ct = ',Ct:2:4,' pf   Cf = ',Cf:2:4,' pf   Cp = ',Cp:2:4,' pf') ;
writeln('New Ct Rod Length = ',l/2.54:2:3,' in   New gap = ',gap/2.54:2:3,' in   Rod + gap = ',(l+gap)/2.54:2:3,' in')  ;
writeln ;

{Estimate resonator Q and filter loss}

resonatorQ := rodQ(rprime,b,l,lamda0,F) ;

writeln ('Estimated resonator Q = ',resonatorQ:4:1,'    Filter Q = ',0.5*Qf:5:1) ;

SuMg := 0.0 ;							{sum prototype g values}
for j := 1 to N do
	SuMg := SuMg + g[j] ;

Loss := 4.343 * (0.5*Qf/resonatorQ) * SuMg ;	{estimated filter loss}

writeln('Loss = ',Loss:3:2,' dB') ;
writeln ;

{It is possible to simulate an array of coupled rods such as an interdigitated
filer with pairs of single coupled lines.  QUCS has a single pair coupled model. 
The newer QucsStudio has a neat tuning function using sliders.  Ct can
be "tuned."

For details of this simulation technique, see: 

Denig, Carl, "Using Microwave CAD Programs to Anallyze Microstrip 
Interdigital Filters," Microwave Journal, March 1989, pp 147-152.}

Zs := 2.0 * Zoslab ;

writeln('Data for simulation model using coupled lines') ;
writeln ;
writeln('Zs = ',Zs:3:2) ;
writeln ;

for j := 0 to N do
begin
	Zoeprime := 1.0 / (1.0 / Zoe[j] - 1.0 / Zs) ;
	Zooprime := 1.0 / (1.0 / Zoo[j] - 1.0 / Zs) ;
	writeln('Zooprime = ',Zooprime:3:2,'  Zoeprime = ',Zoeprime:3:2) ;
end ;

writeln ;

writeln('Ct = ',Ct:2:4,' pf  Bar length = ',10.0*l:2:3,' mm') ;	

readln ;
	
End .


