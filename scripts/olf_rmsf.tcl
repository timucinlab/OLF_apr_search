# Set the path for output file
set output_path "directory_to_output"
set output [open $output_path w]
puts $output "idx resid rmsf_all rmsf_last500"

# Define number of frames in the trajectory
set n [molinfo top get numframes]

# Atom selection for reference and comparison
set reference [atomselect top "protein and segname PROA and name CA" frame 0]
set compare [atomselect top "protein and segname PROA and name CA"]

# Fit the frames to the reference
for {set i 0} {$i < $n} {incr i} {
    $compare frame $i
    set trans_mat [measure fit $compare $reference]
    $compare move $trans_mat
}

# Select protein atoms
set prot [atomselect top "protein and segname PROA and name CA"]

# Define last frame and first frame for the last 5000 frames
set ls [expr {$n-1}]
if {$n > 5000} {
    set fs [expr {$ls-5000}]
} else {
    set fs 0
}

# Calculate RMSF over all frames and the last 5000 frames
set rmsfa [measure rmsf $prot first 0 last $ls step 1]
set rmsfa1 [measure rmsf $prot first $fs last $ls step 1]

# Get residue IDs and write results to the file, including both rmsfa and rmsfa1
set res [$prot get resid]
for {set i 0} {$i < [$prot num]} {incr i} {
    # Write both RMSF values for all frames and the last 5000 frames
    puts $output "[expr {$i+1}] [lindex $res $i] [lindex $rmsfa $i] [lindex $rmsfa1 $i]"
}

# Close the output file
close $output

# Delete all molecules
mol delete all
