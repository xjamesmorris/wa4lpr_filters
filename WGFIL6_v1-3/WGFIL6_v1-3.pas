program                  WaveGuideFilter ;
uses crt, math ;

{created 10/4/88,1.0 working copy, but never released}
{1.1 revision 8/4/89 released as DOS executable at MicrowaveUpdate 89}
{1.2 Edits done to complile with Free Pascal Complier and Geany IDE and
 extended comments for GPL release.  Released on GPL 5/16/2019}
{1.3 Post filter data entry accepts <CR> for current post diameter 
 but you can enter a new post diameters well.  Mimics orginal Turbo 3 
 functionality lost with later compilers. Added 3 dB BW and passband
 return loss calculation for Chebyshev filters. 02/26/2021}


{ Copyright (C) 1989, 2019, 2021  Dennis G. Sweeney  wa4lpr@arrl.net

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
 Versions are available for just about any processor or operating system
 you can think of.  The Ingrated Development Environment used is Geany,
 see www.geany.org.  Geany works with C, Python, Pascal and a number 
 of other programing languages.  It is availible for Linux, Mac OSX, and 
 Winders.  Both FPC and Geany are free under the GPL.
}

type     Datarray  =   Array [0..20] of real ;

var  a                     : real ; {waveguide broadwall dimension }
	 ch					   : char ; {Carriage return test character}
     D                     : real ; {nomralized post diameter D/a  }
     K                 : Datarray ; {nomralized K inverters        }
     X                 : Datarray ; {normalized shunt reactance    }
     Post, Window      : Datarray ; {Post Dia., Iris width         }
     l                 : Datarray ; {resonator cavity length       }
     g                 : Datarray ; {low pass prototype values     }
     phi               : Datarray ; {resonator electrical length   }
     Lg1, Lg2              : real ; {guide wavelenths at f1, f2    }
     Lg0                   : real ; {guide wavelenth at (Lg1+Lg2)/2}
     L0                    : real ; {free space wavelength at fo   }
     fo, f1, f2            : real ; {filter center, lower, and upper freqs}
     wlamda                : real ; {filter frequency mapping parameter}
     s1, s1prime           : real ; {parameters for WG reatance    }
     Xa, Xb                : real ; {shunt, series WG reactance    }
     temp                  : real ; {temperary storage variable    }
     n                  : integer ; {number of filter elements     }
     j                  : integer ; {counter index                 }
     Flag, m            : integer ; {error flag                    }
     BW                    : real ; {bandwidth of filter           }
     Filter             : integer ; {1 = post filter, 2 = iris filter}
     filtype            : integer ; {1 = butter, 2 = chev, 3 = minloss}
     ripple                : real ; {chebychev passband ripple     }
     Another               : char ; {do you wish another design    }
     str1, str2			  :string ; {input strings}

const   cpyrt = 'copyright (C) 1989, 2019, 2021 Dennis G. Sweeney' ;

Function tan(x : real) : real ;  {Tangent function}
begin
  tan := sin(x) / cos(x) ;
end ;

Function Sinh(x:real) : real ;  {Hyperbolic sine function}
begin
  Sinh := 0.5 * (exp(x) - exp(-x)) ;
end ;

Function Arcsinh(x:real) : real ;  {Hyperbolic arcsine function}
begin
   Arcsinh := ln(x + sqrt(sqr(x) + 1.0)) ;
end ;

{LamdaG       returns guide wavelength
              inputs:        f  operating frequency
                             a  guide broad wall dimension }

Function LamdaG(f, a : real) : real ;
begin
  LamdaG := 300.0 / sqrt(sqr(f) - sqr(150.0 / a)) ;
end ;


{SandSprime        returns:   S1, S1prime as global variables
                   imputs :   a  waveguide broadwall dimension
                              L0 free space wavelength at fo  }
                              
{Sec 5-11 of Marcuvitz, Waveguide Handbook}

PROCEDURE SandSprime(a,L0 : real ) ;
var  m                  : integer ;
     aL0, xx            : real    ;
begin
  aL0 := a / L0 ;
  m := 3 ;  S1 := 0.0 ;       {Sec 5-11 of Marcuvitz, Waveguide Handbook}
  repeat                      {Sum the So series}
    xx := 1.0 / sqrt(sqr(m*1.0) - sqr(2.0 * aL0)) - 1/m ;
    S1 := S1 + xx ;
    m := m + 2 ;
  until (xx/S1) < 2.0E-6 ;
  S1 := 2.0 * S1 - 2.0 ;

  m := 3 ; S1prime := 0.0 ;
  repeat                     {Sum the S2 series}
    xx := sqrt(sqr(m * 1.0) - sqr(2.0 * aL0)) - m + 2.0 * sqr(aL0) / m ;
    S1prime := S1prime + xx ;
    m := m + 2 ;
  until (xx/S1prime) < 2.0E-6 ;
  S1prime := -(2.5) +  sqr(L0 / a) * (11.0 / 12.0 - S1prime) ;
end ;

{SinglePostReactance     returns: Xb    series capactive reactance
                                  Xa    shunt inductive reactance
                         inputs:  a     Guide broadwall dimention
                                  D     D/a normalized post diameter
                                  L0    free space wavelength at fo
                                  Lg0   guide wavelength at fo
                global varibles:  S1, S1prime 
                
							Sec 5-11 of Marcuvitz, Waveguide Handbook}
 

PROCEDURE SinglePostReactance(a, L0, Lg0, D : real ; var Xa, Xb : real) ;

var  D2, D4, So, S2     : real ;
     Dpi, aLg0          : real ;
     
begin
  Dpi := D * pi ;  aLg0 := a / Lg0 ;
  D2 := sqr(0.5 * a * Dpi / L0) ;  D4 := sqr(D2) ;
  So := ln(4.0 / Dpi) ;
  S2 := So + S1prime ;
  So := So + S1 ;
  Xb := -aLg0 * sqr(Dpi) / (1.0 + 2.0 * D2 * (S2 + 0.75)) ;
  Xa := 0.5 * (-Xb + aLg0 * (So - D2 - 0.625 * D4 - 2.0 * D4 *
                                     sqr(S2 - 2.0 * So *sqr(L0 / Lg0)))) ;
end ;

{Kinverters produces the K inverters for a shunt reactance waveguide
 filter.  See Matthaei, Young, and Jones, "Microwave Filters....," p 451.
          input                 : wlamda  waveguide freq mapping parameter
                                : n  number of elements in filter
                                : g   array with low pass prototypes
          output                : K   array with inverter values }

PROCEDURE Kinverters(wlamda : real ; n : integer ; g : Datarray ;
                                                        var K : Datarray) ;
                                                        
var  j                  : integer ;

begin
  K[0] := sqrt(0.5 * pi * wlamda / (g[0] * g[1])) ;
  for j := 1 to n - 1 do
    K[j] := 0.5 * pi* wlamda / sqrt(g[j] * g[j+1]) ;
  K[n] := sqrt(0.5 * pi * wlamda / (g[n] * g[n+1])) ;
end ;

{ShuntReactance     input       K: array with K inverter values
                                n: number of element values
                    output      X: array with corresponding shunt
                                   reactances}


PROCEDURE ShuntReactance(n : integer ; K : Datarray ; var X : Datarray) ;

begin
for j := 0 to n do
   X[j] := K[j] / (1 - sqr(K[j])) ;
end ;

{PostDiameter returns:  D    D/a for the given inverter K
                        phi  the electrical length of the K inverter
                input:  a    WG broad wall dimension
                        L0   free space wavelength at fo
                        Lg0  guide wavelength at fo
                        K    value of the desired K inverter}

PROCEDURE PostDiameter(a, L0, Lg0, K : real ; var D, phi : real ; var Flag : integer) ;

var  xnp1, xnm1, Fn : real ;

   Function Fnpost(a, L0, Lg0, K, D : real ; var phi : real) : real ;
   
   var Xa, Xb            : real ;
   
   begin
        SinglePostReactance(a, L0, Lg0, D,Xa,Xb) ;   {Calculate series C and
                                                  shunt L for given D/a}
     phi := -arctan(2.0 * Xa + Xb) - arctan(Xb) ;
     Fnpost := K - Abs(tan(0.5 * phi + arctan(Xb))) ;
   end ;

begin
Flag := 0 ;
  D := 0.125 ; xnm1 := 0.1 ;
  repeat
    Fn := Fnpost(a, L0, Lg0, K, D, phi) ;
    xnp1 := D - Fn * ((D - xnm1) / (Fn - Fnpost(a, L0, Lg0, K, xnm1, phi))) ;
    xnm1 := D ;  D := xnp1 ;
    If d < 0.00001 Then Flag := 1 ;
  until (abs(D - xnm1) < 2.0E-7) or (Flag = 1) ;
  If d > 0.00001 then  xnp1 := Fnpost(a, L0, Lg0, K, D, phi) ;
end ;


{Inductive Window Dimension   a  : broad wall dimenstion of WG
                              L0 : free space wavelength at fo
                              X  : dersired reactance
                     Returns  d  : window width }

Procedure InductiveWindowDem(var d : real; a,L0,X:real) ;

var            xn,xnm1,xnp1    : real ;
               Fn, Fn2         : real ;
               
  Procedure InductiveWindow(var D,B:real) ;
  
  var          d3,d5,s2,s6,s10   : real ;
  
  begin
    d3 := 1 - sqrt(1-sqr(2*a/(3*L0))) ;
    d5 := 1 - sqrt(1-sqr(2*a/(5*L0))) ;
    s2 := sqr(sin(pi*D/2));  s6 := s2*sqr(s2) ;  s10 := s6*sqr(s2) ;
    B := (1/s2) - 1 - (sqr(1-s2)/(1-d3*s6))*
      (3*d3 + 5*d5*(sqr(2 * s2 - 1 + d3 * s6 * (s2 - 2))/
      ((1-d3*s6)*(1-d5*s10)-15*s6*d5*sqr(1-s2))))  ;
  end ;

  Function Fwindow1(var X,D:real) : real ;
  var         B       : real ;
  begin
    InductiveWindow(D,B) ;
    Fwindow1 := (1/X) - B ;
  end ;

  Function Fwindow2(var X,D:real) : real ;
  var         B       : real ;
  begin
    InductiveWindow(D,B) ;
    Fwindow2 := X - 1.0/B ;
  end ;

begin   {procedure InductiveWindow Dimensionm }
		{This proceedure used an interative numberical to y = f(x).  You
		have y which is the reactance and you want x which is the window
		dimension.}
		  
  xn := 0.3 ; xnm1 := 0.35 ;
  Repeat
    If X < 0.5 then
    begin
      Fn  := Fwindow2(X,xn) ;
      Fn2 := Fwindow2(X,xnm1) ;
    end
    else
    begin
      Fn  := Fwindow1(X,xn) ;
      Fn2 := Fwindow1(X,xnm1) ;
    end ;
    xnp1 := xn - Fn*(xn - xnm1) / (Fn - Fn2) ;
    xnm1 := xn ; xn := xnp1 ;
  Until Abs(xn-xnm1) < 2.0E-7 ;
  D := xn
end ;

{Butterworth     returns  g: array of low pass prototype values for
                          Butterworth response
                 input    N: number of elements }

PROCEDURE Butterworth (N : integer ; var g : Datarray) ;

Var       r            :integer ;
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

{Chebyshev       returns  g: array of low pass prototype values for
                          Chebyshev response
                 input    N: number of elements
                          ripple: passband ripple in dB}

PROCEDURE Chebyshev (N : integer; ripple : real ; var g : Datarray) ;

Var        r                     : integer ;
           Theta1, Theta2, nu, e : real ;

Begin
   g[0] := 1.0 ;
   e := sqrt(exp(ripple * ln(10.0) / 10.0) - 1.0) ;
   nu := Sinh(ARCsinh(1/e)/n) ;
   g[1] := 2 * sin(Pi/(2 * N)) / nu ;
   For r := 1 to N-1 Do
   Begin
      Theta1 := (2 * r - 1) * Pi / (2 * N) ;
      Theta2 := (2 * r + 1) * Pi / (2 * N) ;
      g[r + 1] := 4 * sin(Theta1) * sin(Theta2) /
           ((Sqr(nu) + Sqr(sin(r * Pi / N))) * g[r]) ;
  End ;
  IF Odd(N) THEN g[N + 1] := 1  ELSE g[N + 1] := Sqr(e + Sqrt(1 + e*e)) ;
End ;

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

{MinLoss       returns  g: array of low pass prototype values all =1
               input    N: number of elements  }
               
PROCEDURE MinLoss (N : integer; var g : Datarray) ;

var     j                    :integer ;

begin
  for j := 0 to N+1 do
      g[j] := 1.0 ;
end ;

PROCEDURE SignOnMessage ;
begin
  writeln('                Waveguide Filter Systhesis Program  Version 1.3') ;
  writeln('             copyright (C) 1989, 2019, 2021  GPL 3 Dennis G. Sweeney WA4LPR');
  writeln ;  writeln ;
end ;

{Input requests filter type, waveguide broadwall dimension, and
number of elements

       Input:   none
       Returns: n number of elements
                g lowpass prototype values
                a waveguide width
                filter: 1 = post, 2 = iris type filter
                filtype:
                ripple}

Procedure Input(var a : real ; var g : datarray ; var n : integer ;
        var filter : integer; var filtype : integer ; var ripple : real) ;
begin
  repeat
    writeln(' Filter type:   Butterworth   : 1') ;
    writeln('                Chebyshev     : 2') ;
    writeln('                Equal Element : 3') ;
    write('      Enter # of desired type : ') ; readln(filtype) ;
    If filtype = 2 then
    begin
      write('      Enter passband ripple in DB : ') ;
      readln(ripple) ;
    end ;
  until (filtype >= 1) and (filtype <= 3) ;
  writeln ;
  repeat
    writeln(' Filter structure :    Post   : 1') ;
    writeln('                       Iris   : 2') ;
    write('      Enter # of desired type : ') ; readln(filter) ;
  until (filter = 1) or (filter = 2) ;

  writeln ;
  write(' Enter waveguide width (inches) ') ; readln(a) ;
  a := a * 0.0254 ;
  writeln ;

  repeat
  write(' Enter # of elements : ') ; readln(n) ;
  if (n <= 0) or (n >= 21) then
    writeln(' Element number out of range, try again') ;
  until (n > 0) and (n < 21) ;

  Case filtype of
    1: Butterworth(n,g) ;
    2: Chebyshev(n,ripple,g) ;
    3: MinLoss(n,g) ;
  end ;
end ;

BEGIN



SignOnMessage ;             			{Hi there world}
Input(a,g,n,filter,filtype,ripple) ; 	{a = waveguide broad wall
										g = lowpass prototype values
										n = number of elements
										filter = 1 for ports, 2 for iris
										filtype = 1 Butter, 2 Chev, 3 minloss
										ripple = Chebychev ripple}
Repeat
  writeln ;
  repeat
    write(' Center frequency in MHz = ') ;  readln(Fo) ;  {Input center freq}
    temp := 150 / a ;                              {approximate cutoff freq }
    if Fo <= temp then writeln(' Input freqency below WG cut-off, try again')  ;
  until Fo > temp ;
  write(' Bandwidth in MHz = ') ;         readln(BW) ;  {Input bandwidth  }
  
  F1 := fo - bw / 2 ;  F2 := fo + bw / 2 ;

  Lg1 := LamdaG(f1,a)      ;  Lg2    := LamdaG(f2,a) ;
  fo  := sqrt(f1 * f2)     ;  L0     := 300.0 / fo ;
  Lg0 := (Lg1 + Lg2) / 2.0 ;  wlamda := (Lg1 - Lg2) / Lg0 ;

  Kinverters(wlamda,n,g,K) ; {Calculate K inverters}

  SandSprime(a,L0) ;         {Calculate s and s' for shunt susceptance approx}
  If Filter = 1   then       {if Filter = 1 than make a post filter }
  begin
    writeln ;
    writeln('Post Coupled Filter') ; writeln ;
    writeln('Calculated Post dia    Desired post dia') ;

    j := 0 ;     {Calculate post diameter, cavity electrical length}
    repeat
      PostDiameter(a, L0, Lg0, K[j], post[j], phi[j],Flag) ;
      If flag = 1 then
      begin
        writeln('post[',j,'] too small') ;     
      end
      else
      begin
        post[j] := post[j] * a / 0.0254 ;    {post dia in inches}
        write('Post[',j,'] = ',post[j]:7:4,'           ') ;
        
      {Not elegant but it works! Allows <CR> to take current value or
      you can change the value of the post and enter it.  Will print 
      "input error" if not a numerical input}     
        
		repeat
		until Keypressed ;
		ch := readkey ;
		if ch <> #13 then
		begin
			write(ch) ;
			readln(str1) ;        
			str2 := ch + str1 ;    
			val(str2,temp,m) ;			
			if m <> 0 then 
			begin
				writeln('Input error') ;
				temp := post[j] ;
			end ;
			post[j] := temp ;
	    end	
	    else
	       writeln ;	
		
		 	
        D := post[j] * 0.0254 / a ;
        SinglePostReactance(a, L0, Lg0, D,Xa,Xb) ;
        phi[j] := -arctan(2.0 * Xa + Xb) - arctan(Xb) ;
      end ;
      j := j + 1 ;
    until (j = n + 1) or (flag = 1) ;

    for j := 1 to n do     {Calculate cavity physical length in inches}
      l[j] := (0.5 * Lg0 / pi) * (pi + 0.5 * (phi[j-1] + phi[j])) / 0.0254 ;

    ClrScr ;
    
    writeln('Post Coupled Filter    Center: ',fo:9:3,' MHz   BW: ',BW:9:3,' MHz') ;
    write('                       Fractional BW ',((f2-f1)*100.0/fo):6:2,' %') ;
    writeln('   Guide BW ',(wlamda*100.0):6:2,' %') ;
   
    case filtype of
      1: writeln('                       Butterworth response') ;
      2: begin
			writeln('                       Chebychev response, ripple: ',
                                                 ripple:6:2,' DB') ;
		     writeln('                       3 dB BW:  ',(BW*ThreedBBW(N,ripple)):4:2,
				' MHz    Ripple Return Loss: ',RippleRL(ripple):3:1,' dB') ; 
         end ;                                        
      3: writeln('                       Equal element minloss') ;
    end ;
    writeln ;
    If Flag = 1 then
    begin
      writeln('Post coupled filter not realizable.') ;
      writeln('Try a narrower bandwidth.') ;
    end
    else
    begin
      writeln('Post diameter      Cavity length') ;
      write(' ',post[0]:7:4,'"') ;
      If ((post[0] * 0.0254 / a) > 0.25) then
        writeln('                              Caution:  D/a > 0.25')
      else
        writeln ;
      for j := 1 to n do
      begin
        writeln('                  ',l[j]:7:4,'"')  ;
        write(' ',post[j]:7:4) ;
        If ((post[j] * 0.0254 / a) > 0.25) then
          writeln('                              Caution:  D/a > 0.25')
        else
          writeln ;
      end ;
    end ;
    writeln ; write('Do you wish another design (y/n)? ') ; readln(another) ;
    another := upcase(another) ;
  end
  else
  begin              {If postfilter is false make an iris filter}

    shuntreactance(n,K,X) ;   {calculates required shunt reactance}
    for j := 0 to n do
            InductiveWindowDem(window[j],a,L0,X[j]*Lg0/a) ;
                          {calculates iris width from shunt reactance}
    for j := 1 to n do
      l[j] := (Lg0 / (2.0*pi*0.0254)) *
                         (pi - 0.5*(arctan(2.0*X[j-1]) + arctan(2.0*X[j]))) ;
                         {calculates cavity length}
    
    writeln ;
    writeln('Iris Coupled Filter    Center: ',fo:9:3,' MHz    BW: = ',BW:9:3,' MHz') ;
    write('                       Fractional BW ',((f2-f1)*100.0/fo):6:2,' %') ;
    writeln('   Guide BW ',(wlamda*100.0):6:2,' %') ;
    
    case filtype of
      1: writeln('                       Butterworth response') ;
      2: begin
			writeln('                       Chebychev response, ripple: ',
                                                 ripple:6:2,' DB') ;
			writeln('                       3 dB BW:  ',(BW*ThreedBBW(N,ripple)):4:2,
				' MHz    Ripple Return Loss: ',RippleRL(ripple):3:1,' dB') ; 
	     end ;			
                                                 
      3: writeln('                       Equal element minloss') ;
    end ;
    writeln ;
    writeln('Window width      Cavity length') ;
    writeln ;

    writeln('  ',window[0]*(a/0.0254):6:4,'"') ;
    for j := 1 to n do
    begin
      writeln('                  ',l[j]:6:4,'"') ;
      writeln('  ',window[j]*(a/0.0254):6:4,'"') ;
    end ;
    writeln ; write('Do you wish another design (y/n)? ') ; readln(another) ;
    another := upcase(another) ;
  end ;                {End if}
until another = 'N' ;
END.

