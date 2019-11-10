#Principal Component Vectors Analysis

# mv analvec.tcl /home/ariel/progs/vmd-187/vmd/scripts/vmd/
# vi /home/ariel/progs/vmd-187/vmd/scripts/vmd/tclIndex
# add "source /home/ariel/progs/vmd-187/vmd/scripts/vmd/analvec.tcl" at the end of the opened tclIndex file
# type vmd and open VMD Main -> Extensions -> Tk Console and type "analvec" to start

#or simply run a source command en the TkConsole with the path were you have this script and then "analvec"

proc analvec {} {
#   clear
   puts "Principal Component Vectors Analysis\nOUTPUT FILES will look like vec.dat, vec1.tra, vec2.tra...\nType \"analvec_run\" to start\nType \"analvec_show\" to represent the .tra files generated\nFor MORE INFO about the commands type \"analvec_run_info\" or \"analvec_show_info\"\nif you improve this script, please send it back to arielpetruk@gmail.com"

   proc analvec_run_info {} {

      puts "\nTo run the FIRST PROC just type \"analvec_run eigenvector_file vec_number amplitude\"\nanalvec_run can run with any eigenvectors, not just over the backbone eigenvectors\n"

      puts " . \"eigenvector_file\" is the PTRAJ OUTPUT whos first 3 lines look like:\n\n Eigenvector file: COVAR\n 4404 4404\n   44.33531   53.97282   51.87623   43.32636   54.88519   51.33993   43.66719\n\n   and the calculated eigenvectors are separeted samething like this:\n\n   -0.00924    0.00859   -0.00816   -0.00907    0.01049   -0.00972   -0.00693\n    0.01410\n ****\n    2   244.66501\n    0.00813    0.00020    0.00093    0.00998    0.00178   -0.00001    0.01084\n\n   Changes in this file, could generate some problems to run the script.\n"

      puts " . \"vec_number\" is the AMOUNT OF VECTORS that will be calculated.\n   Their eigenvectors and eigenvalues must be in the ptraj output.\n"

      puts " . \"amplitude\" is the EIGENVECTOR AMPLITUDE that will be used to represent the eigenmodes.\n   Try with 20 to start.\n"

   }

   proc analvec_show_info {} {

      puts "\nTo run the SECOND PROC just type \"analvec_show vec_number examined_residues_pdb\".\nThis proc generate the backbone.pdb file and load the .tra files into VMD.\n"

      puts " . \"vec_number\" is the AMOUNT OF EIGENMODES that will be represented. The value must be equal or lower than de amount of eigenmodes calculated in the previous step.\n"

      puts " . \"examined_residues_pdb\" is a PDB FILE with just the analyzed residues.\n   This file may have the complete residues, not just the backbone.\n"

      puts "BE AWARE: if the ptraj out eigenvetors were not just the backbone (N, CA, C and O), it will be necessary to change lines 91 and 106 of this script for the right representation.\n"

   }

   proc analvec_run {eigenvector_file vec_number amplitude} {

      set ptraj_file [exec grep ** $eigenvector_file]
        exec sed "1d" $eigenvector_file > temp
        exec echo " Eigenvector file: COVAR" > temp2
        exec cat  temp2 temp > temp3
        set ptraj_file [exec grep ** temp3]
        set hn [lindex $ptraj_file 3]
	exec rm temp temp2 temp3
      set lptraj_file [llength $ptraj_file]
      #head number = resid number * backbone atom number * 3 (xyz)
      #set hn [lindex $ptraj_file 7]
      #total vec number
      set vectot [expr {($lptraj_file - 2) / ($hn + 3) - 1 }]

      #contribution
      set sum_autov 0
      set lautov {}
      for {set i [expr {$hn + 7}]} {$i < $lptraj_file} {incr i [expr {$hn + 3}]} {
        #add autovectors
        set sum_autov [expr {$sum_autov + [lindex $ptraj_file $i]}]
        #insert each autovector in a list
        set lautov [linsert $lautov end [lindex $ptraj_file $i]]
      }
      set av_number 0
      set sum_av_value 0
      append out_av "NM\t%\tsum%\n"
      foreach av $lautov {
         set av_number [expr {$av_number + 1}]
         set av_value [expr {$av * 100 / $sum_autov}]
         set sum_av_value [expr {$sum_av_value + $av_value}]
         append out_av [format "%.0f\t%.2f\t%.2f\n" $av_number $av_value $sum_av_value]
      }
      #save the autovector contribution in a file
      set bety [open vec_x100.dat w]
      puts $bety $out_av
      close $bety

      for {set i 0} {$i < $vec_number} {incr i} {

         set out {}
         append out "\n"

         for {set j -$amplitude} {$j < [expr {$amplitude+1}]} {incr j 2} {

            set anal {}

            for {set k 5} {$k < [expr {5+$hn}]} {incr k} {

               set cov [lindex $ptraj_file $k]
               set vec [lindex $ptraj_file [expr {$k+($hn+3)*($i+1)}]]
               set anal [linsert $anal end [format "%.3f" [expr {$cov + $vec * $j}]]]

            }

            for {set l 0} {$l < $hn} {incr l 10} {

               for {set m 0} {$m < 10} {incr m} {
                  set a$m [lindex $anal [expr {$l + $m}]]
               }

               append out "$a0 $a1 $a2 $a3 $a4 $a5 $a6 $a7 $a8 $a9\n"

            }

         }

         set bety [open vec[expr {$i+1}].tra w]
         puts $bety $out
         close $bety

         puts "vec[expr {$i+1}].tra end"

      }

       puts "$vec_number vectors were analyzed"
   }

   proc analvec_show {vec_number examined_residues_pdb} {

      #make backbone pdb
      mol load pdb $examined_residues_pdb
      set selbackbone [atomselect top "backbone"]
      animate write pdb backbone.pdb sel $selbackbone top
      mol delete top

      #load vec.tra
      for {set i 1} {$i <= $vec_number} {incr i} {

         mol load pdb backbone.pdb
         animate delete all
         mol addfile vec${i}.tra waitfor all type crd

         if {$i > 1} {mol off top} else {}

         mol delrep 0 top

         mol representation Ribbons
         mol color ColorID $i
         mol selection {all}
         mol material Opaque
         mol addrep top

      }

      axes location off
      animate speed 0.7
      animate style rock
      animate forward

      puts "check your creation"

   }

}

