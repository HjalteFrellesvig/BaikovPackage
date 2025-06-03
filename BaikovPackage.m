(*   This is version 4.2 of BaikovPackage.                                             *)
(*   Last edited on the 3rd of June 2025.                                              *)
(*   For documentation, see arXiv:2412.01804:                                          *)
(*   "The Loop-by-Loop Baikov Representation - Strategies and Implementation"          *)
(*   Published in: JHEP 04 (2025) 111                                                  *)

BeginPackage[ "BaikovPackage`"]

BaikovStandard::usage = "BaikovStandard[] computes the ingredients for a standard Baikov parametrization.";

BaikovLBL::usage = "BaikovLBL[] computes the ingredients for a loop-by-loop Baikov parametrization.";
    
BaikovCombine::usage = "BaikovCombine[comp] computes a Baikov parametrization given ingredients comp.";

MakeExtraDPRules::usage = "MakeExtraDPRules[comp, tau] computes rules for extra dot-products up to tensor-degree tau, given ingredients comp.";

BPKinematics::usage = "BPKinematics[] returns the kinematic variables of the problem.";
  
SetBPprint::usage = "Sets the value of the debug-variable BPprint.";
GetBPprint::usage = "Returns the value of the debug-variable BPprint.";
SetFactorFinal::usage = "Sets the value of the variable BPFactorFinal.";
GetFactorFinal::usage = "Returns the value of the variable BPFactorFinal.";
SetDPresult::usage = "Sets the value of the variable BPDPresult.";
GetDPresult::usage = "Returns the value of the variable BPDPresult.";
SetExpandAndTest::usage = "Sets the value of the variable BPExpandAndTest.";
GetExpandAndTest::usage = "Returns the value of the variable BPExpandAndTest.";

Internal::usage = "Internal is the list of loop-momenta. LBL parametrizes from right to left.";
External::usage = "External is the list of independent momenta external to the whole integral.";
Propagators::usage = "Propagators is the list of (genuine) propagators.";
PropagatorsExtra::usage = "PropagatorsExtra is the list of additional propagator-type objects needed for the Baikov parametrization.";
Replacements::usage = "Replacements is the list of kinematic replacement rules.";
  
x::usage = "x is a variable reserved for Baikov variables.";   
d::usage = "d is a variable reserved for the spacetime dimensionality.";   
DP::usage = "DP is a variable reserved for dotproducts.";
DPx::usage = "DPx is a variable reserved for dotproducts when used as integration variables.";

Begin[ "Private`"]
  
  
  
BPprint = False;
  
SetBPprint = Function[{y},
  If[y===True||y===False,
    BPprint = y;
  ]];
    
GetBPprint = Function[{}, BPprint];
  
  
BPFactorFinal = True;
  
SetFactorFinal = Function[{y},
  If[y===True||y===False,
    BPFactorFinal = y;
  ]];
    
GetFactorFinal = Function[{}, BPFactorFinal];


BPDPresult = False;
  
SetDPresult = Function[{y},
  If[y===True||y===False,
    BPDPresult = y;
  ]];
    
GetDPresult = Function[{}, BPDPresult];


BPExpandAndTest = True;
  
SetExpandAndTest = Function[{y},
  If[y===True||y===False,
    BPExpandAndTest = y;
  ]];
    
GetExpandAndTest = Function[{}, BPExpandAndTest];



PRINT[y___] := If[BPprint, Print[y]];
  
wrongstring = "\nLikely cause is a wrong or incomplete choice of propagators.";


  
MakeDPRules = Function[{ks, ps},
  Module[{dprules1, dprules2},
  dprules1 = {};
  Do[Do[
    AppendTo[dprules1, (ks[[i]]*ks[[j]] -> DP[ks[[i]], ks[[j]]])];
      , {j, i, Length[ks]}], {i, 1, Length[ks]}];
  dprules2 = Flatten[Table[Table[ks[[i]]*ps[[j]] -> DP[ks[[i]], ps[[j]]], {j, Length[ps]}], {i, Length[ks]}]];
  Join[dprules1, dprules2]]];
  

BPKinematics = Function[{}, Complement[ Variables[ Join[Map[#[[2]]&, Replacements], Propagators, PropagatorsExtra] ], Union[Internal, External]]];
    
    
    
    
    
BaikovLBL = Function[{},
  Module[{AllPropagators, PP, RInternal, LL, dprules, bigdplist, whogoeswhere, useddplist, bigjacdet, bigdetlist, bigpowlist, prefaclist, bigxrules, bigelist, result, props, jlist, kk, qq, gg, ktoq, EE, curdplist, propsexpanded, linear1, zpos, gsol, zpos2, qtok, linear2, xlhs, xeqs, toxrules, solveddps, linear3, curfacpows, curfac, cc1, cc2, cc3, cc4, jacmat1, jacmat2, grammat, baikmat, linear4, gramdet, baikdet, jacdet1, jacdet2, gp1, gp2, bigxrulesfinal, bfrfrhs, endlabel, fromxrules, DPxlist, qqofk, DPxrhs, DPxdpl, DPxsol, DPxret, ExpandAndTestLoc},
  
  ExpandAndTestLoc = BPExpandAndTest;

  result = 0;
  LL = Length[Internal];
  RInternal = Reverse[Internal];
  AllPropagators = Join[Propagators, PropagatorsExtra];
  PP = Length[AllPropagators];
  dprules = MakeDPRules[Union[Internal], External];
  bigdplist = Map[#[[2]] &, dprules];
  PRINT["dprules = ", dprules];
  PRINT["bigdplist = ", bigdplist];
    
  whogoeswhere = Table[0, {PP}];
  Do[Do[
    If[Not[FreeQ[AllPropagators[[i]], RInternal[[j]]]],
      whogoeswhere[[i]] = j;
      Break[]];
    , {j, LL}], {i, PP}];
  PRINT["whogoeswhere = ", whogoeswhere];
    
  useddplist = {};
  bigdetlist = {};
  bigpowlist = {};
  prefaclist = {};
  bigjacdet = 1;
  bigxrules = {};
  bigelist = {};
  DPxlist = {};
    
  (*loop-loop*)
  Do[
    kk = RInternal[[i]];
    PRINT["\n  ", kk];
     
    props = {};
    jlist = {};
    Do[
      If[whogoeswhere[[j]] == i,
        AppendTo[props, AllPropagators[[j]]];
        AppendTo[jlist, j]];
    , {j, PP}];
    PRINT["props = ", props];
    EE = Length[props] - 1;
     
    curdplist = {};
    Do[
      If[Not[FreeQ[bigdplist[[j]], kk]], AppendTo[curdplist, bigdplist[[j]]]];
    , {j, Length[bigdplist]}];
    curdplist = Complement[curdplist, useddplist];
    PRINT["curdplist = ", curdplist];
     
    xlhs = Expand[Expand[props] /. dprules /. Replacements];
    xeqs = Table[xlhs[[i]] == x[jlist[[i]]], {i, Length[xlhs]}];
    PRINT["xeqs = ", xeqs];
    toxrules = Quiet[Solve[xeqs, curdplist]];
    If[toxrules === {},
      Print["xeqs = ", xeqs];
      Print["curdplist = ", curdplist];
      Print["ERROR! x-system could not be solved!", wrongstring];
      Goto[endlabel];
    ];
    toxrules = toxrules[[1]];
    toxrules = Union[toxrules];
    PRINT["toxsol = ", toxrules];
    If[Length[toxrules] != EE + 1,
      Print["ERROR! Discrepency in x-solution!", wrongstring];
      Goto[endlabel];
    ];
    solveddps = Map[(#[[1]]) &, toxrules];
    PRINT["solveddps = ", solveddps];
     
    ktoq = {kk -> qq};
    propsexpanded = Expand[Expand[props /. ktoq] /. Replacements];
     
    linear1 = Table[Expand[Coefficient[propsexpanded[[i]], qq, 1]], {i, Length[propsexpanded]}];
    PRINT["linear1 = ", linear1];
     
    zpos = Position[linear1, 0, {1}];
    If[Length[zpos] > 1,
      Print["oh no! The momenta seem linearly dependent.", wrongstring, "\n linear1 = ", linear1];
      Goto[endlabel];
    ];
    If[Length[zpos] == 0,
      ktoq = {kk -> qq + gg};
      propsexpanded = Expand[Expand[props /. ktoq] /. Replacements];
      linear1 = Table[Expand[Coefficient[propsexpanded[[i]], qq, 1]], {i, Length[propsexpanded]}];
      PRINT["linear1 = ", linear1];
      quadratic1 = Table[Expand[Coefficient[propsexpanded[[i]], qq, 2]], {i, Length[propsexpanded]}];
      PRINT["quadratic1 = ", quadratic1];
      fq = 0;
      Do[
        If[quadratic1[[j]]=!=0,
          fq=j;
          Break[];
        ];
        ,{j,Length[quadratic1]}];
      PRINT["fq = ", fq];
      If[fq==0,
        Print["There seem to not be any quadratic props in ", kk, ". We can't handle it.", wrongstring, "\n quadratic1 = ",quadratic1];
        Goto[endlabel];
      ];
      gsol = Solve[linear1[[fq]] == 0, gg][[1]];
      ktoq = ktoq /. gsol;
      propsexpanded = Expand[Expand[propsexpanded /. gsol] /. Replacements];
      linear1 = Expand[Expand[linear1 /. gsol] /. Replacements];
      PRINT["linear1 = ", linear1];
      zpos2 = Position[linear1, 0, {1}];
      If[Length[zpos2] != 1,
        Print["OH NO! The momenta seem linearly dependent.", wrongstring, "\n linear1 = ", linear1];
        Goto[endlabel];
      ];
    ];
    qtok = Solve[kk == (kk /. ktoq), qq][[1]];
    linear2 = DeleteCases[linear1, 0];
    PRINT["linear2 = ", linear2];
     
    linear3 = 0*linear2;
    Do[
      curfac = 1;
      curfacpows = FactorList[cc1*cc2*linear2[[j]]];
      Do[
        If[Intersection[Variables[cc3*cc4*curfacpows[[k, 1]]], Join[Internal, External]] =!= {},
          curfac *= (curfacpows[[k, 1]]^curfacpows[[k, 2]]);
        ];
      , {k, Length[curfacpows]}];
      linear3[[j]] = curfac;
    , {j, Length[linear3]}];
    PRINT["linear3 = ", linear3];
     
    grammat = Expand[KroneckerProduct[linear3, linear3]] /. dprules /. Replacements;
    
    qqofk = qq/.qtok;
    
    AppendTo[DPxlist, DPx[qqofk, qqofk]];
    Do[
      AppendTo[DPxlist, DPx[qqofk, linear3[[j]]]];
      ,{j,Length[linear3]}];
     
    linear4 = Prepend[linear3, qqofk];
    baikmat = Expand[Expand[KroneckerProduct[linear4, linear4]] /. dprules /. Replacements];
     
    jacmat1 = Table[Table[D[xlhs[[i]], solveddps[[j]]], {j, Length[solveddps]}], {i, Length[xlhs]}];
    jacmat2 = Table[Table[D[baikmat[[1, i]], solveddps[[j]]], {j, Length[solveddps]}], {i, Length[baikmat[[1]]]}];
    
    baikmat = baikmat /. toxrules;
      
    grammat = Factor[grammat];
    baikmat = Factor[baikmat];
    (* NEW!!!!!!!! *)
     
    PRINT["grammat = ", grammat];
    PRINT["baikmat = ", baikmat];
    PRINT["jacmat1 = ", jacmat1];
    PRINT["jacmat2 = ", jacmat2];
     
    If[grammat === {},
      gramdet = 1,
      gramdet = Det[grammat]];
    baikdet = Det[baikmat];
    
    If[ExpandAndTestLoc,
      gramdet = Expand[gramdet];
      baikdet = Expand[baikdet];
    ];
    
    bigdetlist = Join[{gramdet, baikdet}, bigdetlist];
     
    bigdetlist = bigdetlist //. toxrules;
   
    If[ExpandAndTestLoc,
      bigdetlist = Map[Expand, bigdetlist];
      If[Apply[And, Map[FreeQ[bigdetlist, #] &, curdplist]],
        PRINT["all is OK for ", kk],
      
        Print["all is NOT OK for ", kk, ".", wrongstring];
        Print["curdplist = ", curdplist];
        Print["bigdetlist = ", bigdetlist];
        Print["toxrules = ", toxrules];
        Print["props = ", props, "   ", "curext = ", linear3];
        Goto[endlabel];
      ];
    ];
     
    gp1 = (EE - d + 1)/2;
    gp2 = (d - EE - 2)/2;
    bigpowlist = Join[{gp1, gp2}, bigpowlist];
    
    If[BPDPresult,
      jacdet1 = 1;
      jacdet2 = 1;
      ,
      jacdet1 = Det[jacmat1];
      jacdet2 = Det[jacmat2];
    ];
    PRINT["jacdets = ", jacdet1, "  ", jacdet2];
     
    prefaclist = PrependTo[prefaclist, (-I*Pi^(-EE/2)/Gamma[(d - EE)/2])];
    bigjacdet = bigjacdet*jacdet2/jacdet1;
     
    bigxrules = Join[bigxrules, toxrules];
    PrependTo[bigelist, linear3];
     
    useddplist = Join[useddplist, curdplist];
  , {i, Length[RInternal]}];
  
  If[Length[DPxlist] != Length[AllPropagators],
    Print["DPxlist has the wrong length. This is too strange.\nDPxlist = ", DPxlist];
    Goto[endlabel];
  ]
    
  PRINT["\n  finally"];
  
  PRINT["DPxlist = ", DPxlist];
  
  If[BPDPresult,
    DPxrhs = Expand[Expand[DPxlist/.{DPx:>Times}]/.dprules/.Replacements];
    DPxdpl = Union[Cases[DPxrhs, DP[___], All]];
    DPxsol = Quiet[Solve[Table[DPxlist[[i]] == DPxrhs[[i]], {i,Length[DPxlist]}], DPxdpl]];
    If[DPxsol==={},
      Print["There seem to be no solution of DPxsol! \nDPxlist = ",DPxlist];
      Goto[endlabel];
    ];
    DPxsol = DPxsol[[1]];
    fromxrules = Table[x[i] -> Expand[Expand[AllPropagators[[i]]]/.dprules/.Replacements/.DPxsol], {i,Length[AllPropagators]}];
    bigdetlist = bigdetlist/.fromxrules;
  ];
  
  If[BPDPresult,
    DPxret = {DPxlist, fromxrules};
    ,
    DPxret = {DPxlist};
  ];
    
  If[ExpandAndTestLoc,
    bigdetlist = Map[Expand, bigdetlist];
  ];
  If[BPFactorFinal, bigdetlist = Factor[bigdetlist]];
  PRINT["bigdetlist = ", bigdetlist];
    
  AppendTo[prefaclist, bigjacdet];

  bfrfrhs = Map[#[[2]] &, bigxrules];
  bfrfrhs = Expand[bfrfrhs //. bigxrules];
  If[BPDPresult, bfrfrhs = Expand[bfrfrhs //. fromxrules]];
  bigxrulesfinal = Table[bigxrules[[i, 1]] -> bfrfrhs[[i]], {i, Length[bfrfrhs]}];

  PRINT["bigxrulesfinal = ", bigxrulesfinal];
    
  result = {bigdetlist, bigpowlist, prefaclist, {dprules, bigxrulesfinal, bigelist, DPxret}};
  Label[endlabel];
  result]];


  



BaikovStandard = Function[{},
  Module[{RInternal, LL, EE, DPrules, AllPropagators, PP, PPexp, AllPropagatorsExpanded, DPRHS, ijacdet, ToBaik, curE, gram1, gp1, curexte, gram2, gp2, remdps, bigdetlist, bigpowlist, prefaclist, Elists, result, endlabel, fromxrules, DPxret},

  result = 0;
  RInternal = Reverse[Internal];
  LL = Length[RInternal];
  EE = Length[External];
  DPrules = MakeDPRules[Internal, External];
  AllPropagators = Join[Propagators, PropagatorsExtra];
  PP = Length[AllPropagators];
  PPexp = (LL*(LL + 1)/2 + LL*EE);
  If[PP != PPexp,
    Print["Wrong number of propagators for standard Baikov!\n L(L+1)/2+LE = ", PPexp, " ,  P = ", PP, wrongstring];
    Goto[endlabel];
  ];
    
  AllPropagatorsExpanded = Expand[Expand[AllPropagators] //. DPrules /. Replacements];
  DPRHS = Union[Cases[AllPropagatorsExpanded, DP[___], All]];
  If[Length[DPRHS] != PP,
    Print["Unmatched dotproducts:", AllPropagatorsExpanded, wrongstring];
    Goto[endlabel];
  ];
  ijacdet = Det[Table[Table[Coefficient[AllPropagatorsExpanded[[i]], DPRHS[[j]]], {j, PP}], {i, PP}]];
  If[ijacdet === 0,
    Print["Extremely bad: ijacdet = 0", wrongstring];
    Goto[endlabel];
  ];
  ToBaik = Solve[Table[AllPropagatorsExpanded[[i]] == x[i], {i, Length[AllPropagatorsExpanded]}], DPRHS][[1]];
    
  curE = Length[External];
  If[curE == 0,
    gram1 = 1, 
    gram1 = Det[KroneckerProduct[External, External] /. DPrules /. Replacements /. ToBaik];
  ];
  gp1 = (curE - d + 1)/2;
  curexte = Join[Internal, External];
  gram2 = Det[KroneckerProduct[curexte, curexte] /. DPrules /. Replacements /. ToBaik];
  gp2 = (d - curE - LL - 1)/2;
  remdps = Union[Cases[{gram1, gram2}, DP[___], All]];
  If[remdps =!= {},
    Print["Oh No! The ", RInternal[[i]], "-loop has unmatched dot products: ", remdps, wrongstring];
    Goto[endlabel];
  ];
  bigdetlist = {gram1, gram2};
  bigpowlist = {gp1, gp2};
  If[BPDPresult,
    ijacdet = 1;
  ];
  prefaclist = {((-I)^LL*Pi^((LL - PPexp)/2)/Product[Gamma[(d + 1 - EE - l)/2], {l, 1, LL}]), 1/ijacdet};
  Elists = {External};
  
  If[BPDPresult,
    fromxrules = Table[x[i] -> AllPropagatorsExpanded[[i]]/.{DP:>DPx} ,{i,PP}];
    bigdetlist = bigdetlist/.fromxrules;
    ToBaik = ToBaik/.fromxrules;
  ];
  
  DPxret = {DPRHS/.{DP:>DPx}};
  If[BPDPresult, AppendTo[DPxret, fromxrules]];
  
  If[BPFactorFinal, bigdetlist = Factor[bigdetlist]];
    
  result = {bigdetlist, bigpowlist, prefaclist, {DPrules, ToBaik, Reverse[Elists], DPxret}};
  Label[endlabel];
  result]];

  
  
  

BaikovCombine = Function[{xx}, Module[{res},
  res = Apply[Times, xx[[3]]];
  Do[
    res *= xx[[1, i]]^xx[[2, i]];
    , {i, Length[xx[[1]]]}];
  res]];
  
  
  
  

MakeExtraDPRules = Function[{inm, nin},
  Module[{resultrules, subresult, dpendlabel, murules, rul1, rul2, alldps, useddps, RInternal, kk, newcandi, pe1, pi1, pi2, allnewcandi, eta, temprules, tensorbasis, pebasis, pebasislist, peb, cfl, ansatz, alist, rhs, lhs,
      eqs, rfl, minusrfl, curminusrfl, rfrules, rfrulesextra, eqssol, eqssolorig, resrhs, reslhs, lhsp, lhsk, rfrulesrev, rfrulesrevextra, ccou, currhs, resultrulespre (*,P,g,mu,a,BR,cc*) },
      
    resultrules = {};
    If[Not[IntegerQ[nin]] || nin < 0,
      Print["N has to be a positive integer"];
      Goto[dpendlabel]];
    If[nin >= 8,
      Print["Not implemented for N>=8"];
      Goto[dpendlabel]];
    If[nin == 0,
      Goto[dpendlabel]];
    
    subresult = Table[{}, {nin}];
    
    murules = {P[x1_, x3_]*P[x2_, x3_] :> BR[x1*x2], P[x1_, x3_]^2 :> BR[x1^2], P[x1_, x2_]*g[x3___, x2_, x4___] :> P[x1, x3, x4], 
      g[x1___, x0_, x2___]*g[x3___, x0_, x4___] :> g[x1, x2, x3, x4], g[x1___, x0_, x2___]^2 :> g[x1, x2, x1, x2], g[x1_, x1_] :> d};
    rul1 = inm[[4, 1]];
    rul2 = inm[[4, 2]];
    alldps = Table[rul1[[i, 2]], {i, Length[rul1]}];
    useddps = Table[rul2[[i, 1]], {i, Length[rul2]}];
    RInternal = Reverse[Internal];
    allnewcandi = {};
    temprules = {};
    
    (*This is the k-loop*)
    Do[
      kk = RInternal[[i]];
      PRINT["k-loop: ", kk, " (", i, "/", Length[RInternal], ")"];
      newcandi = Join[Table[kk*RInternal[[j]], {j, i, Length[RInternal]}], Table[kk*External[[j]], {j, 1, Length[External]}]] /. rul1;
      newcandi = Complement[newcandi, useddps];
      allnewcandi = Join[allnewcandi, newcandi];
      If[Length[newcandi] > 0,
        pe1 = newcandi /. {DP[kk, x1_] :> x1, DP[x1_, kk] :> x1};
        pi1 = inm[[4, 3, Length[Internal] + 1 - i]];
        pi2 = Append[Map[P[#, mu[0]] &, pi1], P[eta, mu[0]]] //. {P[x1_ + x2_, x0_] :> P[x1, x0] + P[x2, x0]};
      
        (*this is the n-loop*)
        Do[
          PRINT["  n-loop: ", n, "/", nin];
          tensorbasis = {1};
          pebasis = {1};
          Do[
            tensorbasis = Flatten[KroneckerProduct[tensorbasis, (pi2 /. {mu[0] -> mu[u]})]];
            pebasis = Union[Flatten[KroneckerProduct[pebasis, pe1]]];
            , {u, 1, n}];
          pebasislist = {};
          Do[
            peb = {};
            cfl = FactorList[pebasis[[j]]];
            Do[
              Do[
                AppendTo[peb, cfl[[u, 1]]];
              , {iii, cfl[[u, 2]]}];
            , {u, 2, Length[cfl]}];
            AppendTo[pebasislist, peb];
          , {j, Length[pebasis]}];
          If[n >= 6,
            tensorbasis = tensorbasis /. {P[eta, x1_]*P[eta, x2_]*P[eta, x3_]*P[eta, x4_]*P[eta, x5_]*P[eta, x6_] :> {g[x1, x6]*g[x2, x5]*g[x3, x4], g[x1, x5]*g[x2, x6]*g[x3, x4], g[x1, x6]*g[x2, x4]*g[x3, x5],
                g[x1, x4]*g[x2, x6]*g[x3, x5], g[x1, x5]*g[x2, x4]*g[x3, x6], g[x1, x4]*g[x2, x5]*g[x3, x6], g[x1, x6]*g[x2, x3]*g[x4, x5], g[x1, x3]*g[x2, x6]*g[x4, x5], g[x1, x2]*g[x3, x6]*g[x4, x5], 
                g[x1, x5]*g[x2, x3]*g[x4, x6], g[x1, x3]*g[x2, x5]*g[x4, x6], g[x1, x2]*g[x3, x5]*g[x4, x6], g[x1, x4]*g[x2, x3]*g[x5, x6], g[x1, x3]*g[x2, x4]*g[x5, x6], g[x1, x2]*g[x3, x4]*g[x5, x6]}};
            tensorbasis = Flatten[tensorbasis];
          ];
          If[n >= 4,
            tensorbasis = tensorbasis /. {P[eta, x1_]*P[eta, x2_]*P[eta, x3_]*P[eta, x4_] :> {g[x1, x2]*g[x3, x4], g[x1, x3]*g[x2, x4], g[x1, x4]*g[x2, x3]}};
            tensorbasis = Flatten[tensorbasis];
          ];
          tensorbasis = tensorbasis /. {P[eta, x1_]*P[eta, x2_] :> g[x1, x2]};
          tensorbasis = DeleteCases[tensorbasis /. {P[eta, x1_] -> 0}, 0];
       
          If[Length[tensorbasis] != 0,
            ansatz = Sum[a[i]*tensorbasis[[i]], {i, Length[tensorbasis]}];
            alist = Union[Cases[ansatz, a[_], All]];
            rhs = Expand[ansatz*tensorbasis] //. murules /. Replacements /. rul1 /. rul2 /. {BR[x1_] :> x1};
            rhs = Map[Collect[#, a[_], Factor] &, rhs];
            rfl = {};
            Do[Do[
              cfl = FactorList[Coefficient[rhs[[ii]], alist[[jj]]]];
              Do[
                If[cfl[[kk, 1, 0]] == Plus, AppendTo[rfl, cfl[[kk, 1]]]];
              , {kk, Length[cfl]}];
            , {jj, Length[alist]}], {ii, Length[rhs]}];
            rfl = Union[rfl];
            minusrfl = Expand[-rfl];         (* TO HERE *)
        
            rfrulesrev = {};
            rfrules = {};
            ccou = 0;
            Do[
              curminusrfl = Take[minusrfl, ii - 1];
              If[FreeQ[curminusrfl, rfl[[ii]]],
                ++ccou;
                AppendTo[rfrulesrev, (cc[ccou] -> rfl[[ii]])];
                rfrules = Flatten[Join[rfrules, {x1_*rfl[[ii]]^x2_ -> x1*cc[ccou]^x2, x1_*rfl[[ii]] -> x1*cc[ccou], x1_*minusrfl[[ii]]^x2_ -> x1*(-1)^x2*cc[ccou]^x2, x1_*minusrfl[[ii]] -> -x1*cc[ccou]}]];
              ];
            , {ii, Length[rfl]}];
        
            PRINT["    ", rfrulesrev];
        
            rhs = rhs //. rfrules;
            PRINT["    ", "rhs:\n", rhs];
          ];
       
          Do[
            PRINT["    p-loop ", pebasislist[[u]], " (", u, "/", Length[pebasislist], ")"];
            lhsp = Apply[Times, Table[P[pebasislist[[u, j]], mu[j]], {j, Length[pebasislist[[u]]]}]];
            lhsk = Apply[Times, Table[P[kk, mu[j]], {j, n}]];
            reslhs = Expand[lhsk*lhsp] //. murules /. Replacements /. rul1 /. rul2 /. {BR[x1_] :> x1};
        
            If[Length[tensorbasis] === 0,
              resrhs = 0;
              PRINT["      No tensor basis -> no solving"];
         
              ,
         
              lhs = Expand[lhsp*tensorbasis] //. murules /. Replacements /. rul1 /. rul2 /. {BR[x1_] :> x1};
              lhs = Map[Factor, lhs];
              lhs = lhs //. rfrules;
         
              rfl = {};
              Do[
                cfl = FactorList[lhs[[ii]]];
                Do[
                  If[cfl[[kk, 1, 0]] == Plus, AppendTo[rfl, cfl[[kk, 1]]]];
                , {kk, Length[cfl]}];
              , {ii, Length[lhs]}];
              rfl = Union[rfl];
              minusrfl = Expand[-rfl];
         
              rfrulesextra = {};
              rfrulesrevextra = {};
              ccou = Length[rfrulesrev];
              Do[
                curminusrfl = Take[minusrfl, ii - 1];
                If[FreeQ[curminusrfl, rfl[[ii]]],
                  ++ccou;
                  AppendTo[rfrulesrevextra, (cc[ccou] -> rfl[[ii]])];
           
                  rfrulesextra = Flatten[Join[rfrulesextra, {x1_*rfl[[ii]]^x2_ -> x1*cc[ccou]^x2, x1_*rfl[[ii]] -> x1*cc[ccou], x1_*minusrfl[[ii]]^x2_ -> x1*(-1)^x2*cc[ccou]^x2, x1_*minusrfl[[ii]] -> -x1*cc[ccou]}]];
                ];
              , {ii, Length[rfl]}];
         
              lhs = lhs //. rfrulesextra;
         
              rfrules = Join[rfrules, rfrulesextra];
              rfrulesrev = Join[rfrulesrev, rfrulesrevextra];
         
              PRINT["      ", rfrulesrev];
              PRINT["      ", "lhs:\n", lhs];
         
              eqs = Table[lhs[[j]] == rhs[[j]], {j, Length[rhs]}];
              PRINT["      ", "Solving (size=", Length[rhs], ")"];
              eqssolorig = Solve[eqs, alist];
              PRINT["      ", "Solving Done"];
              If[Length[eqssolorig] != 1,
                Print["Something is wrong with the a-solution: Length[sols] = ", Length[eqssolorig]];
              ];
              PRINT["      ", "bytecount = ", ByteCount[eqssolorig]];
              PRINT["      ", "Factoring sol"];
              eqssol = Table[0, {Length[alist]}];
              
              Do[
                PRINT["        ", j, "/", Length[alist], "   ", ByteCount[eqssolorig[[1, j, 2]]]];
                currhs = Factor[eqssolorig[[1, j, 2]]];
                eqssol[[j]] = (eqssolorig[[1, j, 1]] -> currhs);
              , {j, Length[alist]}];
              eqssolorig = 0;
              PRINT["      ", "Factoring sol done"];
              PRINT["      ", "bytecount = ", ByteCount[eqssol]];
              PRINT["      ", "      eqssol = ", eqssol];
         
              resrhs = Expand[lhsk*(ansatz /. eqssol)] //. murules /. Replacements /. rul1 /. rul2 /. {BR[x1_] :> x1};
              PRINT["      ", "Factoring"];
              resrhs = Factor[resrhs];
              PRINT["      ", "Factoring Done"];
              resrhs = resrhs /. rfrulesrev;
              (* resrhs=Factor[resrhs]; *)
            ];
            AppendTo[subresult[[n]], (reslhs -> resrhs)];
          , {u, Length[pebasislist]}];
        , {n, nin, 1, -1}];
      ];
    , {i, Length[RInternal]}];
    resultrulespre = Flatten[Reverse[subresult]];
    PRINT["The last steps"];
    resultrules = Table[0, {i, Length[resultrulespre]}];
    Do[
      currhs = resultrulespre[[i, 2]];
      While[Not[FreeQ[currhs, DP[___]]],
        PRINT["  ", i, "/", Length[resultrulespre]];
        currhs = Expand[currhs] /. resultrulespre;
        (*currhs=Factor[currhs];*)
      ];
      resultrules[[i]] = (resultrulespre[[i, 1]] -> currhs);
    , {i, Length[resultrulespre]}];
    PRINT["Done!"];
    Label[dpendlabel];
    resultrules
  ]];

  

End[]

EndPackage[]
