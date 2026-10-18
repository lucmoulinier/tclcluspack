# Appels repetes dans un meme interpreteur (comme dans Ordali) :
# verifie la stabilite et la croissance memoire.
# Usage : tclsh8.6 stress.tcl <builddir> ?nbIterations?
set dir [file dirname [info script]]
set argv [list [lindex $argv 0] __aucun__ [lindex $argv 1]]
source [file join $dir cases.tcl]

proc rss {} {
    set f [open /proc/self/status]
    set s [read $f]
    close $f
    regexp {VmRSS:\s+(\d+)} $s -> kb
    return $kb
}

set n [lindex $argv 2]
if {$n eq ""} {set n 200}
set appels {
    {Tclcluspack $A -dt alignment -cm bionj -nbc secator}
    {Tclcluspack $A -dt alignment -cm kmeans -nbc dpc -fd 0.001}
    {Tclcluspack $A -dt alignment -cm hierar -nbc secator -fd 0.001}
    {Tclcluspack $C -dt coordinates -cm kmeans -nbc dpc -wc}
    {Tclcluspack $C -dt coordinates -cm hierar -nbc secator -wc}
    {Tclcluspack $C -dt coordinates -cm mm -nbc bic}
}
set rss0 [rss]
for {set i 1} {$i <= $n} {incr i} {
    foreach cmd $appels {
        if {[catch {eval $cmd} res]} {
            puts "iteration $i : ERREUR $res ($cmd)"
            exit 1
        }
    }
    if {$i == 10} {set rss10 [rss]}
}
puts "OK $n iterations x [llength $appels] appels ; RSS debut $rss0 Ko, apres 10 iter. $rss10 Ko, fin [rss] Ko"
