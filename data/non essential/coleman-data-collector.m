// Note: this records the *full* output of the ColemanData intrinsic, resulting in hundreds of MB of data.
// Use short-coleman-data-collector.m instead unless you absolutely need to store the value of one of the big fields.
AttachSpec("coleman.spec");
load "data/quartic-test-data.m";

coleman_data_format:=recformat<Q,p,N,g,W0,Winf,r,Delta,s,G0,Ginf,e0,einf,delta,basis,quo_map,
  integrals,F,f0list,finflist,fendlist,Nmax,red_list_fin,red_list_inf,minpolys,cpm,subspace,
  ordinary,frobmatb0r>;

coleman_data_list := [* *];
output_file := "data/coleman-test-data.m";
if not OpenTest(output_file, "r") then
  fprintf output_file, "coleman_data_format := %m;\n\ncoleman_data_list := [* *];\n\n", coleman_data_format;
  fprintf output_file, "_<x> := PolynomialRing(RationalField());\n_<z> := LaurentSeriesRing(PolynomialRing(RationalField()));\n\n";
end if;

for i := 1 to #test_data_list do
  test := test_data_list[i];
  printf "Computing Coleman data for Q = %o...", test`Q;
  data := ColemanData(test`Q, test`p, test`N : useU:=true, basis0:=test`basis0, basis1:=test`basis1, basis2:=test`basis2);
  Append(~coleman_data_list, data);
  out := Sprintf("coleman_data_list[%o] := %m;\n\n", i, data);
  out := ReplaceString(out, "\n ! ", "\n "); // hack to handle bug in Magma string formatting
  Write(output_file, out);
  printf " done.\n";
end for;
