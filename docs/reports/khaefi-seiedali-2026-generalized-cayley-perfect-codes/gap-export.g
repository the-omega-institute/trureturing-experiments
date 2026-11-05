LoadPackage("smallgrp");
ExportData:=function()
 local n,idx,G,els,auts,invol,subs,aut,s,H,mask,row,col;
 for n in [1..16] do
  for idx in [1..NumberSmallGroups(n)] do
   G:=SmallGroup(n,idx); els:=Elements(G);
   Print("GROUP ",n," ",idx,"\n");
   for row in [1..n] do
    for col in [1..n] do Print(Position(els,els[row]*els[col])," "); od; Print("\n");
   od;
   auts:=Elements(AutomorphismGroup(G));
   invol:=Filtered(auts,aut->ForAll(els,g->Image(aut,Image(aut,g))=g));
   for aut in invol do
    Print("AUT"); for g in els do Print(" ",Position(els,Image(aut,g))); od; Print("\n");
   od;
   subs:=AllSubgroups(G);
   for H in subs do
    if Size(H)<n then
     mask:=0; for g in Elements(H) do mask:=mask+2^(Position(els,g)-1); od;
     Print("SUB ",mask,"\n");
    fi;
   od;
   Print("END\n");
  od;
 od;
end;
ExportData();
QUIT;
