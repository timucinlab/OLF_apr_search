# Set the path for output file
set output_path "path_to_output.dat"
set output [open $output_path w]
puts $output "idx resid rgyr"

# Define number of frames in the trajectory
set n [molinfo top get numframes]

# Atom selection for reference and comparison
set reference [atomselect top "protein and name CA" frame 0]
set compare [atomselect top "protein and name CA"]

# Fit the frames to the reference
for {set i 0} {$i < $n} {incr i} {
    $compare frame $i
    set trans_mat [measure fit $compare $reference]
    $compare move $trans_mat
}

# Select protein atoms
set prot [atomselect top "protein and name CA and resid 381 326 380 429 428 477 478 328 476 271 329 269"]


# calculation loop
for {set i 0} {$i < $n} {incr i} {
  $prot frame $i
  set rgyr [measure rgyr $prot]
  puts $output "$i $rgyr"
}


# Close the output file
close $output

# Delete all molecules
mol delete all
