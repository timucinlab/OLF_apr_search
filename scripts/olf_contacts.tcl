set n [molinfo top get numframes]
# Open the output file
set output [open "ure6.contacts.dat" w]
puts $output "frame bladeAB_atomcontact bladeBC_atomcontact bladeCD_atomcontact bladeDE_atomcontact bladeEA_atomcontact bladeAB_rescontact bladeBC_rescontact bladeCD_rescontact bladeDE_rescontact bladeEA_rescontact "

# Define selections
set sel0 [atomselect top "(protein and resid 268 to 317) and within 4 of (protein and resid 328 to 359)"]
set sel1 [atomselect top "(protein and resid 328 to 359) and within 4 of (protein and resid 380 to 422)"]
set sel2 [atomselect top "(protein and resid 380 to 422) and within 4 of (protein and resid 428 to 467)"]
set sel3 [atomselect top "(protein and resid 428 to 467) and within 4 of (protein and resid 474 to 501 255 to 259)"]
set sel4 [atomselect top "(protein and resid 474 to 501 255 to 259) and within 4 of (protein and resid 268 to 317)"]

# Get the number of frames in the trajectory
set n [molinfo top get numframes]

# Loop through each frame
for {set i 0} {$i < $n} {incr i} {
    $sel0 frame $i
    $sel0 update
    set sel00 [lsort -unique [$sel0 get {resid resname}]]

    $sel1 frame $i
    $sel1 update
    set sel11 [lsort -unique [$sel1 get {resid resname}]]

    $sel2 frame $i
    $sel2 update
    set sel22 [lsort -unique [$sel2 get {resid resname}]]

    $sel3 frame $i
    $sel3 update
    set sel33 [lsort -unique [$sel3 get {resid resname}]]

    $sel4 frame $i
    $sel4 update
    set sel44 [lsort -unique [$sel4 get {resid resname}]]

    # Write results to output file
    puts $output "$i [$sel0 num] [llength $sel00] [$sel1 num] [llength $sel11] [$sel2 num] [llength $sel22] [$sel3 num] [llength $sel33] [$sel4 num] [llength $sel44]"
}

# Close the output file
close $output
