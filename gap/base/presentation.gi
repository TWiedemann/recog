#############################################################################
##
##  This file is part of recog, a package for the GAP computer algebra system
##  which provides a collection of methods for the constructive recognition
##  of groups.
##
##  This files's authors include Torben Wiedemann.
##
##  Copyright of recog belongs to its developers whose names are too numerous
##  to list here. Please refer to the COPYRIGHT file for details.
##
##  SPDX-License-Identifier: GPL-3.0-or-later
##
##
##  This file provides code for standard presentations of some groups
##  (which are not covered by the classicpres package).
##
#############################################################################

# Returns the presentation of S_n from
# [CM80] "Generators and relations for discrete groups", (6.21).
# The same presentation appears in
# [BLN+03] "A black-box group algorithm for recognizing finite symmetric 
# and alternating groups, I", (2.1),
# except that j <= n/2 in [BLN+03] whereas j <= n-2 in [CM80]. It is not
RECOG.SnPresentation := function(n)
    local F, rels;
    if n=1 then
        return FreeGroup(0);
    elif n=2 then
        F := FreeGroup(1);
        return F / [ F.1^2 ];
    else
        F := FreeGroup(2);
        rels := [ F.1^2, F.2^n, (F.1*F.2)^(n-1), (F.1*(F.1^F.2))^3 ];
        rels := Concatenation(
            rels, List([2..n-2], j -> (F.1*(F.1^(F.2^j)))^2)
        );
        return F / rels;
    fi;
end;

#
RECOG.AnPresentation := function(n)
    local F, rels;
    F := FreeGroup(2);
    rels := [ F.2^(n-2), F.1^3 ];
    if IsOddInt(n) then
        Add(rels, (F.2*F.1)^n);
        rels := Concatenation(
            rels,
            List([1..QuoInt(n-3, 2)], k -> (F.1*F.2^(-k)*F.1*F.2^k)^2)
        );
    else
        Add(rels, (F.2*F.1)^(n-1));
        rels := Concatenation(
            rels,
            List([1..QuoInt(n-2, 2)], k -> (F.1^((-1)^k)*F.2^-k*F.1*F.2^k)^2)
        );
    fi;
    return F / rels;
end;
