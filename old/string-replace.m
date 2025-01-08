freeze;

intrinsic ReplaceCharacter(s::MonStgElt,c::MonStgElt,d::MonStgElt) -> MonStgElt
{ Efficiently replace every occurrence of the character c in s with the string d (c must be a single character but d need not be). }
    require #c eq 1: "The second parameter must be a single character (string of length 1).";
    t := Split(s,c:IncludeEmpty);
    if s[#s] eq c then Append(~t,""); end if; // add empty string at the end which Split omits
    return Join(t,d);
end intrinsic;

intrinsic ReplaceString(s::MonStgElt,c::MonStgElt,d::MonStgElt) -> MonStgElt
{ Greedily replace each occurrence of string c in s with the string d. }
    require #c ge 1: "The string to be replaced cannot be empty.";
    m := #c;
    t := "";
    n := Index(s,c);
    while n gt 0 do
        t cat:= s[1..n-1] cat d;
        s := s[n+m..#s];
        n := Index(s,c);
    end while;
    return t cat s;
end intrinsic;

