# Cas de test pour tclcluspack.
# Usage : tclsh8.6 cases.tcl <builddir> ?nomDuCas?
lappend auto_path [lindex $argv 0]
package require tclcluspack
expr {srand(42)}

# 3 nuages gaussiens 2D bien separes
proc nuages {{n 20}} {
    set L {}
    set k 0
    foreach {cx cy} {0 0 10 10 20 0} {
        for {set i 0} {$i < $n} {incr i} {
            lappend L [list p[incr k] [expr {$cx + rand()*2}] [expr {$cy + rand()*2}]]
        }
    }
    return $L
}

# alignement : 3 familles de sequences (liste nom seq nom seq ...)
proc alignement {{n 6}} {
    set L {}
    set k 0
    set aa ACDEFGHIKLMNPQRSTVWY
    foreach base {ACDEFGHIKLMNPQRSTVWYACDEFGHIKL MNPQRSTVWYACDEFGHIKLMNPQRSTVWY WYVTSRQPNMLKIHGFEDCAYWVTSRQPNM} {
        for {set i 0} {$i < $n} {incr i} {
            set s $base
            for {set m 0} {$m < 3} {incr m} {
                set p [expr {int(rand()*[string length $s])}]
                set s [string replace $s $p $p [string index $aa [expr {int(rand()*20)}]]]
            }
            lappend L seq[incr k] $s
        }
    }
    return $L
}

# matrice de distances euclidiennes {nom d1 ... dn}
proc distances {L} {
    set D {}
    foreach a $L {
        set row [list [lindex $a 0]]
        foreach b $L {
            lappend row [expr {hypot([lindex $a 1]-[lindex $b 1], [lindex $a 2]-[lindex $b 2])}]
        }
        lappend D $row
    }
    return $D
}

# 3 gaussiennes en 1 dimension (moyennes 0, 3.5, 7 ; ecart-type ~1) :
# l'ancienne penalite BIC (p.ln(n) au lieu de p.ln(n)/2) ne trouvait qu'un groupe
proc gauss1d {{n 50}} {
    set L {}
    set k 0
    foreach m {0 3.5 7} {
        for {set i 0} {$i < $n} {incr i} {
            # somme de 12 uniformes - 6 : approximativement N(0,1)
            set x -6.0
            for {set u 0} {$u < 12} {incr u} {set x [expr {$x + rand()}]}
            lappend L [list g[incr k] [expr {$m + $x}]]
        }
    }
    return $L
}

# compositions en acides amines (20 frequences de somme 1, covariance
# singuliere) : 3 familles de 60 sequences de composition differente
proc compositions {{n 60}} {
    set L {}
    set k 0
    foreach famille {0 1 2} {
        for {set i 0} {$i < $n} {incr i} {
            set v {}
            for {set a 0} {$a < 20} {incr a} {
                set base [expr {($a % 3 == $famille) ? 3.0 : 1.0}]
                lappend v [expr {$base * (0.8 + 0.4*rand())}]
            }
            set s [tcl::mathop::+ {*}$v]
            lappend L [list c[incr k] {*}[lmap x $v {expr {$x/$s}}]]
        }
    }
    return $L
}

set C [nuages]
set G [gauss1d]
set P [compositions]
set P15 [compositions 5]
set A [alignement]
set D [distances [nuages 8]]
set N {((a:0.1,b:0.1):2.0,(c:0.1,d:0.1):2.0,((e:0.1,f:0.1):2.0,fictive_sequence:3.0):0.5);}
# memes nuages, coordonnees x1000 (kkmeans : bandwidth independant de l'echelle)
set C1000 [lmap p $C {list [lindex $p 0] [expr {[lindex $p 1]*1000}] [expr {[lindex $p 2]*1000}]}]
set cases [dict create \
    coord_kmeans_dpc      {Tclcluspack $C -dt coordinates -cm kmeans -nbc dpc} \
    coord_kmeans_dpc_wc   {Tclcluspack $C -dt coordinates -cm kmeans -nbc dpc -wc} \
    coord_kmeans_3        {Tclcluspack $C -dt coordinates -cm kmeans -nbc 3} \
    coord_kmeans_gap      {Tclcluspack $C -dt coordinates -cm kmeans -nbc gap} \
    coord_hierar_secator  {Tclcluspack $C -dt coordinates -cm hierar -nbc secator -wc} \
    coord_hierar_3        {Tclcluspack $C -dt coordinates -cm hierar -nbc 3} \
    coord_hierar_dpc      {Tclcluspack $C -dt coordinates -cm hierar -nbc dpc} \
    coord_mm_aic          {Tclcluspack $C -dt coordinates -cm mm -nbc aic} \
    coord_mm_bic          {Tclcluspack $C -dt coordinates -cm mm -nbc bic} \
    coord_mm_3            {Tclcluspack $C -dt coordinates -cm mm -nbc 3} \
    gauss1d_mm_bic        {Tclcluspack $G -dt coordinates -cm mm -nbc bic} \
    gauss1d_mm_aic        {Tclcluspack $G -dt coordinates -cm mm -nbc aic} \
    compo_mm_bic          {Tclcluspack $P -dt coordinates -cm mm -nbc bic} \
    compo_mm_aic          {Tclcluspack $P -dt coordinates -cm mm -nbc aic} \
    fewind_mm_bic         {Tclcluspack $P15 -dt coordinates -cm mm -nbc bic} \
    fewind_mm_aic         {Tclcluspack $P15 -dt coordinates -cm mm -nbc aic} \
    err_fewind_mm_diag0   {Tclcluspack $P15 -dt coordinates -cm mm -nbc bic -diag 0} \
    coord_emm_3           {Tclcluspack $C -dt coordinates -cm emm -nbc 3} \
    coord_kkmeans_3       {Tclcluspack $C -dt coordinates -cm kkmeans -nbc 3} \
    coord_kkmeans_3_x1000 {Tclcluspack $C1000 -dt coordinates -cm kkmeans -nbc 3} \
    ali_kkmeans_3         {Tclcluspack $A -dt alignment -cm kkmeans -nbc 3} \
    err_kkmeans_toomany   {Tclcluspack $C -dt coordinates -cm kkmeans -nbc 61} \
    coord_kmm_3           {Tclcluspack $C -dt coordinates -cm kmm -nbc 3} \
    ali_bionj_secator     {Tclcluspack $A -dt alignment -cm bionj -nbc secator} \
    ali_kmeans_dpc_fd     {Tclcluspack $A -dt alignment -cm kmeans -nbc dpc -fd 0.001} \
    ali_hierar_secator_fd {Tclcluspack $A -dt alignment -cm hierar -nbc secator -fd 0.001} \
    ali_bionj_2           {Tclcluspack $A -dt alignment -cm bionj -nbc 2} \
    ali_bionj_secator_fd  {Tclcluspack $A -dt alignment -cm bionj -nbc secator -fd 0.001} \
    ali_kmeans_3          {Tclcluspack $A -dt alignment -cm kmeans -nbc 3} \
    ali_mm_bic            {Tclcluspack $A -dt alignment -cm mm -nbc bic} \
    dist_kmeans_3         {Tclcluspack $D -dt distances -cm kmeans -nbc 3} \
    dist_hierar_secator   {Tclcluspack $D -dt distances -cm hierar -nbc secator -wc} \
    dist_bionj_secator    {Tclcluspack $D -dt distances -cm bionj -nbc secator} \
    dist_kmeans_3_fd      {Tclcluspack $D -dt distances -cm kmeans -nbc 3 -fd 0.5} \
    newick_secator        {Tclcluspack $N -dt newick -nbc secator} \
    err_newick_nofictive  {Tclcluspack {((a:0.1,b:0.1):2.0,(c:0.1,d:0.1):2.0,(e:0.1,f:0.1):2.0);} -dt newick -nbc secator} \
    err_noargs            {Tclcluspack} \
    err_nodt              {Tclcluspack $C -cm kmeans -nbc 3} \
    err_missing_value     {Tclcluspack $C -dt coordinates -cm} \
    err_bad_option        {Tclcluspack $C -dt coordinates -foo} \
    err_bad_number        {Tclcluspack $C -dt coordinates -cm kmeans -nbc abc} \
    err_bad_row           {Tclcluspack {{a 1 2} {b 1} {c 3 4}} -dt coordinates -cm kmeans -nbc 2} \
    err_not_numeric       {Tclcluspack {{a 1 2} {b 1 x} {c 3 4}} -dt coordinates -cm kmeans -nbc 2} \
    err_ali_odd           {Tclcluspack {s1 ACDE s2} -dt alignment -cm bionj -nbc secator} \
    err_mm_dpc            {Tclcluspack $C -dt coordinates -cm mm -nbc dpc} \
    err_nocm              {Tclcluspack $C -dt coordinates -nbc 3} \
    err_bad_dt            {Tclcluspack $C -dt xxx -cm kmeans -nbc 3} \
]

set only [lindex $argv 1]
dict for {name cmd} $cases {
    if {$only ne "" && $only ne $name} continue
    set code [catch {eval $cmd} res]
    if {$code} {
        puts "$name: ERREUR: $res"
    } else {
        set tailles [lmap g $res {llength $g}]
        puts "$name: OK [llength $res] groupes, tailles $tailles"
    }
}
