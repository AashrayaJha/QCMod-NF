load "data/Correspondence_250.m";
load "data/Correspondence_400.m";

bad_indices:=[[ 1, 5 ],[ 1, 6 ],[ 2, 4 ],[ 2, 6 ],[ 3, 4 ],[ 3, 5 ]];

for i in [1..6] do
    for j in [1..6] do
        if not [i,j] in bad_indices then 
            a:=AK_patch_3_250[i,j]-AK_patch_3_400[i,j];
            if not a eq 0 then
                printf "The values don't agree at %o,%o\n",i,j;
            end if;    
        end if;    
    end for;
end for;        
