#
gap> START_TEST("presentation.tst");
gap> for n in [1..7] do
>  if Size(RECOG.SnPresentation(n)) <> Size(SymmetricGroup(n)) then
>    Print("Presentation for symmetric group S_", n, " has incorrect size\n");
>  fi;
>  if Size(RECOG.AnPresentation(n)) <> Size(AlternatingGroup(n)) then
>    Print("Presentation for alternating group A_", n, " has incorrect size\n");
>  fi;
> od;
gap> STOP_TEST("verification.tst");
