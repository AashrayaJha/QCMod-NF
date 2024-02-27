AttachSpec("QCMod.spec");
//load "data/NF-example-test-data.m";
load "data/NF-example-coleman-data_40.m";
data:=data;
//test:=test_NF;
QpPoints:=data`Qppoints_40;
P:=QpPoints[9];
// import "applications.m" :Fp_points,Qp_points;
// //import "singleintegrals.m" :is_very_bad,is_bad;
//import "singleintegrals.m" :tiny_integrals_on_basis,is_very_bad,find_bad_point_in_disk,tadicprec,is_bad,local_coord;
//import "misc.m": alg_approx_Qp, alg_dep_powerseries;
x0:=P`x; b:=P`b; Q:=data`Q; p:=data`p; N:=data`N; basis:=data`basis; r:=data`r; W0:=data`W0; Winf:=data`Winf;
d:=Degree(Q); lc_r:=LeadingCoefficient(r); W:=Winf*W0^(-1); Qp:=Parent(x0); K:=BaseRing(BaseRing(Q)); v:=data`v;
prec:=20;

if is_bad(P,data) and not is_very_bad(P,data) then // on a bad disk P needs to be very bad
    P1:=find_bad_point_in_disk(P,data);  
else
    P1:=P;
end if;
x1:=P1`x;

IPP1,NIPP1:=tiny_integrals_on_basis(P,P1,data:prec:=20);

if prec eq 0 then // no t-adic precision specified
    prec:=tadicprec(data,1);
end if;

Qpt<t>:=LaurentSeriesRing(Qp,prec);
Kt<t>:=LaurentSeriesRing(K,prec);
Zp:=RingOfIntegers(Qp);
Zpt:=LaurentSeriesRing(Zp,prec);

// xt,bt,index:=local_coord(P1,prec,data);

// xtold:=xt;
// btold:=bt;

// //xt:=Kt!alg_dep_powerseries(xt,v);
// btnew:=[Kt|];
// for i:=1 to d do
//   btnew[i]:=Kt!alg_dep_powerseries(btold[i],v);
// end for;
// bt:=Vector(btnew);

// if P1`inf then
//     xt:=1/xt;
//     xt:=Kt!Qpt!xt; 
//     Winv:=W0*Winf^(-1);          
//     bt:=bt*Transpose(Evaluate(Winv,xt));
//     for i:=1 to d do
//       bt[i]:=Kt!(Qpt!bt[i]);
//     end for; 
//   end if;

//   if P1`inf or not is_bad(P1,data) then 
//     denom:=Qpt!(1/Evaluate(r,xt));
//   else
//     Qp:=pAdicField(p,N);
//     Qpx:=PolynomialRing(Qp);
//     rQp:=Qpx!r;
//     zero:=HenselLift(rQp,x1);
//     sQp:=rQp div (Qpx.1-zero);
//     denom:=Kt!Qpt!((Kt!Zpt!(xt-Coefficient(xt,0)))^(-1)*(Kt!Qpt!(1/Evaluate(sQp,xt))));
//   end if;