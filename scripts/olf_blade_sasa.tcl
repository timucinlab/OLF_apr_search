set n [molinfo top get numframes]

set output [open "6M_urea_diSu.sasa.dat" w]
puts $output "frame sasa_protein sasa_bladeA_wrt sasa_bladeA sasa_bladeB_wrt sasa_bladeB sasa_bladeC_wrt sasa_bladeC sasa_bladeD_wrt sasa_bladeD sasa_bladeE_wrt sasa_bladeE"

set sel [atomselect top "protein"]
set sel0 [atomselect top "protein and segname PROA and resid 268 to 317"]
set sel1 [atomselect top "protein and segname PROA and resid 328 to 359"]
set sel2 [atomselect top "protein and segname PROA and resid 380 to 422"]
set sel3 [atomselect top "protein and segname PROA and resid 428 to 467"]
set sel4 [atomselect top "protein and segname PROA and resid 255 to 259 474 to 501"]

# calculation loop
for {set i 0} {$i < $n} {incr i} {
	puts "\t \t progress: $i/$n"
	molinfo top set frame $i
	set sasa [measure sasa 1.4 $sel]
	set sasa0 [measure sasa 1.4 $sel -restrict $sel0]
	set sasa1 [measure sasa 1.4 $sel0]
	set sasa2 [measure sasa 1.4 $sel -restrict $sel1]
	set sasa3 [measure sasa 1.4 $sel1]
	set sasa4 [measure sasa 1.4 $sel -restrict $sel2]
	set sasa5 [measure sasa 1.4 $sel2]
	set sasa6 [measure sasa 1.4 $sel -restrict $sel3]
	set sasa7 [measure sasa 1.4 $sel3]
	set sasa8 [measure sasa 1.4 $sel -restrict $sel4]
	set sasa9 [measure sasa 1.4 $sel4]
	
	

	puts $output "$i $sasa $sasa0 $sasa1 $sasa2 $sasa3 $sasa4 $sasa5 $sasa6 $sasa7 $sasa8 $sasa9"
}
close $output




