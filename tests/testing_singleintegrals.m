AttachSpec("~/GitHub/QCMod-NF/QCMod.spec");
//import "singleintegrals.m": lie_in_same_disk;
//import ""
load "data/NF-example-test-data.m";
data:=test_NF;
Points:=data`Qppoints;
for i in [1..12] do
    P1:=Points[i];
    for j in [i+1..13] do
        P2:=Points[j];
        if lie_in_same_disk(P1,P2,data) then
            print i;
            print j;
        end if;   
    end for;
end for;
print "tested points if they lie in the same residue disk"; // No points lie in the same disk

P:=Points[0]
x0:=P`x; b:=P`b; Q:=data`Q; p:=data`p; N:=data`N; basis:=data`basis; r:=data`r; W0:=data`W0; Winf:=data`Winf;
d:=Degree(Q); lc_r:=LeadingCoefficient(r); W:=Winf*W0^(-1); Qp:=Parent(x0); K:=BaseRing(BaseRing(Q));

if is_bad(P,data) and not is_very_bad(P,data) then // on a bad disk P needs to be very bad
    P1:=find_bad_point_in_disk(P,data);  
else
    P1:=P;
end if;
x1:=P1`x;
