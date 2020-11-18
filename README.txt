(* this are number when i get only yields *)
{\[Gamma]sfit -> 0.204392, \[Mu]Sfit -> 
  0.117832, \[Mu]I3fit -> -0.0138863, \[Mu]Bfit -> 0.772769, 
 Tfit -> 0.0476622, rmaxfit -> 17.501}
{hh -> 0.007462941669361537`};
r = 17.596313428173186`;

(* this are number when i get yields + all errors *)
{\[Gamma]sfit -> 0.24117850465190532`, \[Mu]Sfit -> 
  0.11603383764876621`, \[Mu]I3fit -> -0.013637416049298717`, \
\[Mu]Bfit -> 0.7706168899154414`, Tfit -> 0.046942706032170475`, 
 rmaxfit -> 
  18.375091500381096`}; h0 = 0.007067; r = 18.54387235541713`;

(* this are number when i get yields - all errors *)
{\[Gamma]sfit -> 0.167534, \[Mu]Sfit -> 
  0.120055, \[Mu]I3fit -> -0.0141832, \[Mu]Bfit -> 0.775427, 
 Tfit -> 0.0485537, rmaxfit -> 16.5266}
 
 {hh -> 0.007916736138354174`}; r = 16.558427791559588`;
 
 
 
 dN\[Pi]p = 0.64; N\[Pi]p = 9.32(*-dN\[Pi]p*); 
 dN\[Pi]m = 1.1; N\[Pi]m = 17.1 (*- dN\[Pi]m*);
  dNKp = 0.00897; NKp = 0.0598 (*- dNKp*); 
 dNKm = 8.55 10^-5; NKm = 0.00057 (*- dNKm*);
(*Np=77.46;    dNp=2.;*)
(*Np=77.46+ 28.7;    dNp=2. +0.79;*)
(* \
tutaj blad powinien byc dodany do siebie kwadratow i wyciagniety \
pierwiastek 2.41*)
(* koniec koncow bledy sa dodane do siebie liniowo \
4.13 *)
  dNp = 2. + 0.79 + 1.05 + 0.29; Np = 
 77.46 + 28.7 + 8.65 + 4.55 (*- dNp*);  
(*Np=81.5+30.6;    dNp=12.225+4.6;*)
dN\[CapitalLambda] = 0.013155; N\
\[CapitalLambda] = 0.0877 