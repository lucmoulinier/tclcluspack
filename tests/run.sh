#!/bin/bash
# Lance chaque cas de tests/cases.tcl dans un processus separe (un crash n'arrete pas les autres).
# Usage : tests/run.sh ?builddir? ?cas ...?
d=$(cd "$(dirname "$0")" && pwd)
b=${1:-build-debug}; shift
tclsh=${TCLSH:-/home/moumou/tcl/bin/tclsh8.6}
cases=${*:-$(sed -n 's/^    \([a-z_0-9]*\) *{Tclcluspack.*/\1/p' "$d/cases.tcl")}
for c in $cases; do
    timeout 60 "$tclsh" "$d/cases.tcl" "$b" "$c" > /tmp/cluspack_test.$$ 2>&1
    rc=$?
    out=$(grep -v "^nb clusters\|^meilleurNombre\|^fin du modele" /tmp/cluspack_test.$$ | tail -1)
    printf '%-24s rc=%-3s %s\n' "$c" "$rc" "${out:0:150}"
done
rm -f /tmp/cluspack_test.$$
