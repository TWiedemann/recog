#############################################################################
##
##  This file is part of recog, a package for the GAP computer algebra system
##  which provides a collection of methods for the constructive recognition
##  of groups.
##
##  This file's authors include Torben Wiedemann.
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

# Returns the presentation Fp (i.e. an FpGroup) of S_n from
# [BLN+03] "A black-box group algorithm for recognizing finite symmetric 
# and alternating groups, I", (2.1)
# or an "obvious" presentation of S_n for small n.
# For n=1, Fp is the trivial group.
# For n>1, Fp has a generator Fp.1 which corresponds to (1,2).
# This is the only generator if n=2.
# For n>2, Fp has exactly two generators and Fp.2 corresponds to (1,...,n).
# Note: The claim in [BLN+03] that this presentation can be found in
# [CM80] "Generators and relations for discrete groups"
# is apparently slightly imprecise. In [CM80, (6.26)],
# the same presentation minus the relations
# s^n = 1 and (t*s)^(n-1) but with the additional weaker relation
# s^n = (t*s)^(n-1) appears. Since this presentation defines S_n by [CM80], we
# can add additional relations that are satisfied in S_n without changing the
# presented group, so the presentation in [BLN+03, (2.1)] is correct.
RECOG.SnPresentation := function(n)
    local F, rels, s, t, j;
    # Edge cases for small n
    if n<1 then
        return fail;
    elif n=1 then
        return FreeGroup(0);
    elif n=2 then
        F := FreeGroup(1);
        return F / [ F.1^2 ];
    fi;
    # Generic case
    F := FreeGroup(2);
    t := F.1; # Notation for s and t as in [BLN+03]
    s := F.2;
    rels := [ s^n, t^2, (s*t)^(n-1) ];
    for j in [2..QuoInt(n,2)] do
        Add(rels, (t*s^-j*t*s^j)^2);
    od;
    return F / rels;
end;

# Returns the presentation Fp (i.e. an FpGroup) of A_n from
# [BLN+03] "A black-box group algorithm for recognizing finite symmetric 
# and alternating groups, I", (2.2), (2.3).
# or an "obvious" presentation of A_n for small n.
# For n in [1,2], Fp is the trivial group.
# For n=3, Fp has a single generator Fp.1.
# For n>3, Fp has exactly two generators: Fp.1 corresponds to (1,2,3) and
# Fp.2 corresponds to
# - (3,...,n) if n is odd,
# - (1,2)(3,...,n) if n is even.
# For odd n, this is the same presentation as in
# [CM80] "Generators and relations for discrete groups", (2.2).
# For even n, it is similar to [CM80, (2.3)].
RECOG.AnPresentation := function(n)
    local F, rels, s, t, k;
    # Edge cases for small n
    if n<1 then
        return fail;
    elif n=1 or n=2 then
        return FreeGroup(0);
    elif n=3 then
        F := FreeGroup(1);
        return F / [F.1^3];
    fi;
    # Generic case
    F := FreeGroup(2);
    t := F.1; # Notation for s and t as in [BLN+03]
    s := F.2;
    rels := [ s^(n-2), t^3 ];
    if IsOddInt(n) then
        Add(rels, (s*t)^n);
        for k in [1..QuoInt(n-3, 2)] do
            Add(rels, (t*s^(-k)*t*s^k)^2);
        od;
    else
        Add(rels, (s*t)^(n-1));
        Add(rels, (t^-1*s^-1*t*s)^2);
    fi;
    return F / rels;
end;
