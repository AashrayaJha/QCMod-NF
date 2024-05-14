AttachSpec("coleman.spec");
load "data/quartic-test-data.m";

short_coleman_data_format:=recformat<Q,p,N,g,W0,Winf,r,Delta,s,G0,Ginf,e0,einf,delta,basis,
  integrals,F,Nmax,minpolys,cpm,subspace,ordinary>;

short_coleman_data_list := [* *];
output_file := "data/short-coleman-test-data.m";
if not OpenTest(output_file, "r") then
  fprintf output_file, "short_coleman_data_format := %m;\n\nshort_coleman_data_list := [* *];\n\n", short_coleman_data_format;
  fprintf output_file, "_<x> := PolynomialRing(RationalField());\n\n";
end if;

for i := 1 to #test_data_list do
  test := test_data_list[i];
  printf "Computing Coleman data for Q = %o...", test`Q;
  full_data := ColemanData(test`Q, test`p, test`N : useU:=true, basis0:=test`basis0, basis1:=test`basis1, basis2:=test`basis2);
  data := rec< short_coleman_data_format | >;
  for s in Names(short_coleman_data_format) do
    if assigned full_data``s then
        data``s := full_data``s;
    end if;
  end for;
  Append(~short_coleman_data_list, data);
  out := Sprintf("short_coleman_data_list[%o] := %m;\n\n", i, data);
  out := ReplaceString(out, "\n ! ", "\n "); // hack to handle bug in Magma string formatting
  Write(output_file, out);
  printf " done.\n";
end for;
