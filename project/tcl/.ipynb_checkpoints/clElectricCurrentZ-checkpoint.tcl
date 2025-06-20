# Calculate the current for a trajectory. 
# Results are in "step current(A)"

# Parameters
## dt in s
set dt 5e-11
set startFrame 0
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Output
set outFile Desktop/clcurrent.dat

# System size in A
set lx 68.476
set ly 68.476
set lz 136.401

# Open the output file
set out [open $outFile w]

# Summation loop
set sel [atomselect top "name CLA and z > 50 and z < 70"]
molinfo top set frame $startFrame
set z0 [$sel get z]
$sel set charge -1.60217733e-19
set q [$sel get charge]

# Begin summation
for {set f 1} {$f < $nFrames} {incr f} {
     molinfo top set frame $f
     $sel update
     
    if {[$sel num] == 0} {
        puts $out "$f 0"
        puts "No selected atoms, writing zero current."
        continue
     } 

     $sel set charge -1.60217733e-19
     set q [$sel get charge]
     
     # Get the position data for the current frame
     set z1 [$sel get z]
     
     # Find the displacements in the z-direction
     set dz {}
     if {[llength $z0] == [llength $z1]} {
        foreach r0 $z0 r1 $z1 {
            # Compensate for jumps across the periodic cell
            set z [expr $r1 - $r0]
            if {[expr $z > 0.5 * $lz]} {set z [expr $z - $lz]}
            if {[expr $z < -0.5 * $lz]} {set z [expr $z + $lz]}
            lappend dz $z
            }
    
     } else {
    # Set all dz values based on the assumption that z0 = 70
        foreach r1 $z1 {
            lappend dz 0
        }
    }
     # Compute the average charge*velocity between the two frames in C A/s
     set qvsum [expr [vecdot $dz $q] / $dt]
          
     # Scale by the system size to obtain the z-current in C/s
     set currentZ [expr $qvsum/$lz]
               
     # Write the current
     puts $out "$f $currentZ"
     puts -nonewline [format "FRAME %i: " $f]
     puts "$f $currentZ"
     
     $sel update
    
     # Store the postion data for the next computation.
     set z0 [$sel get z]
}
close $out
