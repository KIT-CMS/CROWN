for MASS_H3 in 280 #240 280 320 360 400 450 500 550 600 700 800 900 1000
do

#for MASS_H2 in 60 70 80 90 100 120 150 170 190 250 300 350 400 450 500 550 600 650 700 800 900 1000 1100 1200 1300 1400 1600 1800 2000 2200 2400 2600 2800
for MASS_H2 in 60 70 75 80 85 90 95 100 110 120 130 150 170 190 250 300 350 400 450 500 550 600 650 700 750 800 850
do
SUM=$((MASS_H2+125))

if [ "$SUM" -gt "$MASS_H3" ]
then
   echo "Skipping "$MASS_H2" "$MASS_H3
   continue      # Skip rest of this particular loop iteration.
fi

mv M${MASS_H3}_h1_M125_tautau_h2_M${MASS_H2}_bb_miniaod_batch2.txt M${MASS_H3}_h1_M125_tautau_h2_M${MASS_H2}_bb_miniAOD.txt


done

done
