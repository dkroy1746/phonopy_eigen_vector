#!/bin/bash

echo -e " I need 1.band.yaml containing eigen vectors and 2. POSCAR-unitcell used in phonon calculation"
echo -e "primitive cell from band.yaml of phonopy" >> primitiv.vasp
echo -e "1.0" >> primitiv.vasp
grep -w -A3 "lattice:" band.yaml |awk 'NR>1{print "     ",$3,$4,$5}'|sed 's/,/ /g' >> primitiv.vasp
awk 'NR==6{print $0}' POSCAR-unitcell >> primitiv.vasp
awk 'NR==7{for(i=1;i<=NF;i++){printf "   %s  ", $i/2}}' POSCAR-unitcell >> primitiv.vasp
echo -e "\nDirect" >> primitiv.vasp
grep -w "coordinates:" band.yaml |awk '{print "     ",$3,$4,$5}'|sed 's/,/ /g' >> primitiv.vasp
a1=$(grep -w -A3 "lattice:" band.yaml |awk 'NR==2{print "     ",$3}'|sed 's/,/ /g')
a2=$(grep -w -A3 "lattice:" band.yaml |awk 'NR==2{print "     ",$4}'|sed 's/,/ /g')
a3=$(grep -w -A3 "lattice:" band.yaml |awk 'NR==2{print "     ",$5}'|sed 's/,/ /g')
b1=$(grep -w -A3 "lattice:" band.yaml |awk 'NR==3{print "     ",$3}'|sed 's/,/ /g')
b2=$(grep -w -A3 "lattice:" band.yaml |awk 'NR==3{print "     ",$4}'|sed 's/,/ /g')
b3=$(grep -w -A3 "lattice:" band.yaml |awk 'NR==3{print "     ",$5}'|sed 's/,/ /g')
c1=$(grep -w -A3 "lattice:" band.yaml |awk 'NR==4{print "     ",$3}'|sed 's/,/ /g')
c2=$(grep -w -A3 "lattice:" band.yaml |awk 'NR==4{print "     ",$4}'|sed 's/,/ /g')
c3=$(grep -w -A3 "lattice:" band.yaml |awk 'NR==4{print "     ",$5}'|sed 's/,/ /g')

i_alpha=$(echo "($b1*$c1+$b2*$c2+$b3*$c3)/(sqrt(($b1)^2+($b2)^2+($b3)^2)*sqrt(($c1)^2+($c2)^2+($c3)^2))" | bc -l)
echo $a1 $a2 $a3 $b1 $b2 $b3 $c1 $c2 $c3
if [ "${i_alpha}" == 0  ]
          then
          alpha=90.0
      else
          alpha=$(echo "scale=5;180+a(sqrt(1-($i_alpha)^2)/$i_alpha)*180/3.14159265359" | bc -l)
fi


i_beta=$(echo "($a1*$c1+$a2*$c2+$a3*$c3)/(sqrt(($a1)^2+($a2)^2+($a3)^2)*sqrt(($c1)^2+($c2)^2+($c3)^2))" | bc -l)
if [ "${i_beta}" == 0  ]
          then
          beta=90.0
      else
          beta=$(echo "scale=5;180+a(sqrt(1-($i_beta)^2)/$i_beta)*180/3.14159265359" | bc -l)
fi

i_gamma=$(echo "($b1*$a1+$b2*$a2+$b3*$a3)/(sqrt(($b1)^2+($b2)^2+($b3)^2)*sqrt(($a1)^2+($a2)^2+($a3)^2))" | bc -l)
if [ "${i_gamma}" == 0  ]
          then
          gamma=90.0
      else
          gamma=$(echo "scale=5;180+a(sqrt(1-($i_gamma)^2)/$i_gamma)*180/3.14159265359" | bc -l)
fi

echo $alpha $beta $gamma

# Write in POSCAR.mcif using primitiv.vasp
echo -e "\n_chemical_name_common\t '$(awk 'NR==1 {print $0}' primitiv.vasp)' " > POSCAR.mcif
echo -e "_cell_length_a\t\t $(awk 'NR==3 {print  sqrt($1^2+$2^2+$3^2)}' primitiv.vasp | bc -l)" >> POSCAR.mcif
echo -e "_cell_length_b\t\t $(awk 'NR==4 {print  sqrt($1^2+$2^2+$3^2)}' primitiv.vasp | bc -l)" >> POSCAR.mcif
echo -e "_cell_length_c\t\t $(awk 'NR==5 {print  sqrt($1^2+$2^2+$3^2)}' primitiv.vasp | bc -l)" >> POSCAR.mcif
echo -e "_cell_angle_alpha\t\t $alpha">> POSCAR.mcif
echo -e "_cell_angle_beta\t\t $beta">> POSCAR.mcif
echo -e "_cell_angle_gamma\t\t $gamma">> POSCAR.mcif

echo -e "\nloop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z" >> POSCAR.mcif

no_of_species=$(awk ' {if ( NR == 6 )  print NF }' primitiv.vasp)

no_of_atoms=$(awk -v j=0 'NR==7{for(i=1;i<=NF;i++){j=$i+ j}; print j}' primitiv.vasp)

NION=0
for (( i=1; i<=$no_of_species; i++ ))
do
       spec[$i]=$(awk -v j="$i" '{ if ( NR == 6 ) print $j }' primitiv.vasp)
       tot_atm_spec[$i]=$(awk -v j="$i" '{ if ( NR == 7 ) print $j }' primitiv.vasp)

            for (( w=1; w<=${tot_atm_spec[$i]}; w++  ))
            do
                    echo -e "${spec[$i]}$w\t $(awk -v x=$(($NION+$w+8)) 'NR==x {print $0}' primitiv.vasp)" >> POSCAR.mcif

            done

       st_no[$i]=$(($NION+1))
       echo "Species $i : " ${spec[$i]} ${tot_atm_spec[$i]} " Starts from atom no. " ${st_no[$i]}
       NION=$(($NION+${tot_atm_spec[$i]}))

done


tail -$no_of_atoms POSCAR.mcif > temp.txt
dmmy=$(($no_of_atoms*4+1))

echo -e -n '\nloop_
_atom_site_moment.label
_atom_site_moment.crystalaxis_x
_atom_site_moment.crystalaxis_y
_atom_site_moment.crystalaxis_z\n'>>  POSCAR.mcif

read -p " Enter the frequency (only the first q-point eigen vector will be considered corresponding to the frequency) " freq

for (( i=1; i<=$no_of_atoms; i++ ))
 do
     awk -v a=$i 'NR==a{printf "%s ",$1}' temp.txt >> temp2.txt

            for (( j=2; j<=4; j++  ))
            do
                k=$(grep -m1 -A$dmmy "frequency:   $freq" band.yaml |grep -A3 "atom $i" | awk -v w=$j 'NR==w{ print $3}'|sed 's/,//g')
              # echo $k $i
               awk -v a=$i -v b=$k 'NR==a {printf " %f ",b}' temp.txt >> temp2.txt
            done
        echo -n -e "\n" >> temp2.txt
 done

cp POSCAR.mcif temp.txt
cat temp.txt temp2.txt > POSCAR.mcif
rm temp.txt temp2.txt
mv POSCAR.mcif POSCAR_$freq.mcif
