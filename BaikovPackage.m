(*This is version 4.0 of BaikovPackage. Last edited on the 26th of August 2024*)

BeginPackage[ "BaikovPackage`"]

BaikovStandard::usage = "BaikovStandard[] computes the ingredients for a standard Baikov parametrization.";

BaikovLBL::usage = "BaikovLBL[] computes the ingredients for a loop-by-loop Baikov parametrization.";
    
BaikovCombine::usage = "BaikovCombine[x] computes a Baikov parametrization given ingredients x.";
  
SetBPprint::usage = "Sets the value of the debug-variable BPprint.";
GetBPprint::usage = "Returns the value of the debug-variable BPprint.";
SetFactorFinal::usage = "Sets the value of the variable BPFactorFinal.";
GetFactorFinal::usage = "Returns the value of the variable BPFactorFinal.";
SetDPresult::usage = "Sets the value of the variable BPDPresult.";
GetDPresult::usage = "Returns the value of the variable BPDPresult.";

Internal::usage = "Internal is the list of loop-momenta. LBL Integrates from left to right.";
External::usage = "External is the list of independent momenta external to the whole integral.";
Propagators::usage = "Propagators is the list of (genuine) propagators.";
PropagatorsExtra::usage = "PropagatorsExtra is the list of additional propagators needed for the Baikov Parametrization.";
Replacements::usage = "Replacements is the list of kinematic replacement rules.";
  
x::usage = "x is a variable reserved for Baikov variables.";   
d::usage = "d is a variable reserved for the spacetime dimentionality.";   
DP::usage = "DP is a variable reserved for dotproducts.";
DPx::usage = "DPx is a variable reserved for dotproducts when used as integration variables.";

Begin[ "Private`"]
  
  
  
BPprint = False;
  
SetBPprint = Function[{y},
  If[y===True||y===False,
    BPprint = y;
  ]];
    
GetBPprint = Function[{}, BPprint];
  
PRINT[y___] := If[BPprint, Print[y]];
  
wrongstring = "\nLikely cause is a wrong or incomplete choice of propagators.";
  
  
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

  
MakeDPRules = Function[{ks, ps},
  Module[{dprules1, dprules2},
  dprules1 = {};
  Do[Do[
    AppendTo[dprules1, (ks[[i]]*ks[[j]] -> DP[ks[[i]], ks[[j]]])];
      , {j, i, Length[ks]}], {i, 1, Length[ks]}];
  dprules2 = Flatten[Table[Table[ks[[i]]*ps[[j]] -> DP[ks[[i]], ps[[j]]], {j, Length[ps]}], {i, Length[ks]}]];
  Join[dprules1, dprules2]]];
    
    
    
    
    
BaikovLBL = Function[{},
  Module[{AllPropagators, PP, RInternal, LL, dprules, bigdplist, whogoeswhere, useddplist, bigjacdet, bigdetlist, bigpowlist, prefaclist, bigxrules, bigelist, result, props, jlist, kk, qq, gg, ktoq, EE, curdplist, propsexpanded, linear1, zpos, gsol, zpos2, qtok, linear2, xlhs, xeqs, toxrules, solveddps, linear3, curfacpows, curfac, cc1, cc2, cc3, cc4, jacmat1, jacmat2, grammat, baikmat, linear4, gramdet, baikdet, jacdet1, jacdet2, gp1, gp2, bigxrulesfinal, bfrfrhs, endlabel, fromxrules, DPxlist, qqofk, DPxrhs, DPxdpl, DPxsol},

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
    baikmat = Expand[baikmat /. toxrules];
     
    PRINT["grammat = ", grammat];
    PRINT["baikmat = ", baikmat];
    PRINT["jacmat1 = ", jacmat1];
    PRINT["jacmat2 = ", jacmat2];
     
    If[grammat === {},
      gramdet = 1,
      gramdet = Det[grammat]];
    baikdet = Det[baikmat];
    
    gramdet = Expand[gramdet];
    baikdet = Expand[baikdet];
    
    bigdetlist = Join[{gramdet, baikdet}, bigdetlist];
     
    bigdetlist = Map[Expand, bigdetlist //. toxrules];
    If[Apply[And, Map[FreeQ[bigdetlist, #] &, curdplist]],
      PRINT["all is OK for ", kk],
      
      Print["all is NOT OK for ", kk, ".", wrongstring];
      Print["curdplist = ", curdplist];
      Print["bigdetlist = ", bigdetlist];
      Print["toxrules = ", toxrules];
      Print["props = ", props, "   ", "curext = ", linear3];
      Goto[endlabel];
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
    fromxrules = Table[x[i] -> (Expand[AllPropagators[[i]]]/.dprules/.Replacements/.DPxsol), {i,Length[AllPropagators]}];
    bigdetlist = bigdetlist/.fromxrules;
  ];
    
  bigdetlist = Map[Expand, bigdetlist];
  If[BPFactorFinal, bigdetlist = Factor[bigdetlist]];
  PRINT["bigdetlist = ", bigdetlist];
    
  AppendTo[prefaclist, bigjacdet];
  
  If[BPDPresult,
    bigxrulesfinal = fromxrules,
    
    bfrfrhs = Map[#[[2]] &, bigxrules];
    bfrfrhs = Expand[bfrfrhs //. bigxrules];
    bigxrulesfinal = Table[bigxrules[[i, 1]] -> bfrfrhs[[i]], {i, Length[bfrfrhs]}];
  ];
  PRINT["bigxrulesfinal = ", bigxrulesfinal];
    
  result = {bigdetlist, bigpowlist, prefaclist, {dprules, bigxrulesfinal, bigelist, DPxlist}};
  Label[endlabel];
  result]];


  



BaikovStandard = Function[{},
  Module[{RInternal, LL, EE, DPrules, AllPropagators, PP, PPexp, AllPropagatorsExpanded, DPRHS, ijacdet, ToBaik, curE, gram1, gp1, curexte, gram2, gp2, remdps, bigdetlist, bigpowlist, prefaclist, Elists, result, endlabel, fromxrules},

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
    ToBaik = fromxrules;
  ];
    
  If[BPFactorFinal, bigdetlist = Factor[bigdetlist]];
    
  result = {bigdetlist, bigpowlist, prefaclist, {DPrules, ToBaik, Reverse[Elists], DPRHS/.{DP:>DPx}}};
  Label[endlabel];
  result]];

  

BaikovCombine = Function[{xx}, Module[{res},
  res = Apply[Times, xx[[3]]];
  Do[
    res *= xx[[1, i]]^xx[[2, i]];
    , {i, Length[xx[[1]]]}];
  res]];

End[]

EndPackage[]
