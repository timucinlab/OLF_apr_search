# Define atom selections for the two groups
set sel0 [atomselect top "protein and segname PROA and resid 268 to 317"]
set sel1 [atomselect top "protein and segname PROA and resid 328 to 359"]
set sel2 [atomselect top "protein and segname PROA and resid 380 to 422"]
set sel3 [atomselect top "protein and segname PROA and resid 428 to 467"]
set sel4 [atomselect top "protein and segname PROA and resid 474 to 501 255 to 259"]


# Get the total number of frames in the trajectory
set num_frames [molinfo top get numframes]

# Open a file to save the distances
set outfile [open "path_to_output.txt" w]
puts $outfile "Frame\tDistance bladeAB (A)\tDistance bladeBC (A)\tDistance bladeCD (A)\tDistance bladeDE (A)\tDistance bladeEA (A)"

# Loop over each frame to calculate the COM distance
for {set i 0} {$i < $num_frames} {incr i} {
    # Update to the current frame
    animate goto $i
    
    # Calculate the center of mass for each group
    set com0 [measure center $sel0 weight mass]
    set com1 [measure center $sel1 weight mass]
    set com2 [measure center $sel2 weight mass]
    set com3 [measure center $sel3 weight mass]
    set com4 [measure center $sel4 weight mass]
    
    
    # Calculate the distances between each pair of COMs
    set distance_AB [veclength [vecsub $com0 $com1]]
    set distance_BC [veclength [vecsub $com1 $com2]]
    set distance_CD [veclength [vecsub $com2 $com3]]
    set distance_DE [veclength [vecsub $com3 $com4]]
    set distance_EA [veclength [vecsub $com4 $com1]]
    
    # Print and save the frame number and distance
    puts $outfile "$i\t$distance_AB\t$distance_BC\t$distance_CD\t$distance_DE\t$distance_EA"

}

# Close the output file
close $outfile

