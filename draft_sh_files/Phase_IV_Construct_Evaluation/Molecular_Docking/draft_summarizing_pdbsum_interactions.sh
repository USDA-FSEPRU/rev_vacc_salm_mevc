
##Loop for use
for file in 2z7x_TLR1_2_Complex_Any_Model_*_Interactions.txt
do filename=${file%_*}
grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" $file | while read p; do echo $p; done | cut -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > "$filename"_Residues.txt
grep -Fx -f ../2z7x_TLR1_2_Complex_Interactions_Chain_A_Residues.txt "$filename"_Residues.txt > "$filename"_A_Residues_Overlapping.txt
grep -Fx -f ../2z7x_TLR1_2_Complex_Interactions_Chain_B_Residues.txt "$filename"_Residues.txt > "$filename"_B_Residues_Overlapping.txt
Chain_A_Residue_Count="$(wc -l ../2z7x_TLR1_2_Complex_Interactions_Chain_A_Residues.txt | cut -f 1 -d ' ')"
Chain_B_Residue_Count="$(wc -l ../2z7x_TLR1_2_Complex_Interactions_Chain_B_Residues.txt | cut -f 1 -d ' ')"
Query_Chain_A_Count="$(wc -l "$filename"_A_Residues_Overlapping.txt | cut -f 1 -d ' ')"
Query_Chain_B_Count="$(wc -l "$filename"_B_Residues_Overlapping.txt | cut -f 1 -d ' ')"
Chain_A_Percentage=$(echo "scale=2;$Query_Chain_A_Count*100/$Chain_A_Residue_Count" |bc)
Chain_B_Percentage=$(echo "scale=2;$Query_Chain_B_Count*100/$Chain_B_Residue_Count" |bc)
echo $filename $Chain_A_Percentage $Chain_B_Percentage >> 2z7x_TLR1_2_Complex_Any_Loc_Docking_Percentages.txt
done


##Notated Loop
#For every file associated with the TLR5 Any Models
for file in 2z7x_TLR1_2_Complex_Any_Model_*_Interactions.txt

#Extract the filename and save it as a variable
do filename=${file%_*}

#Extract every line that starts with a letter with 0, 1, or 2 spaces in front of it with grep
#Then remove all duplicate and leading white spaces with echo
#Then cut out the Residue Number and Chain Letter for each line witch cut
#Then remove the space in between them with sed
#Finally sort and keep unique values and save them in a file named with filename
grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" $file | while read p; do echo $p; done | cut -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > "$filename"_Residues.txt

#Extract out the Residues in Chain A that match the Reference with grep
grep -Fx -f 2z7x_TLR1_2_Complex_Interactions_Chain_A_Residues.txt "$filename"_Residues.txt > "$filename"_A_Residues_Overlapping.txt

#Extract out the Residues in Chain B that match the Reference with grep
grep -Fx -f 2z7x_TLR1_2_Complex_Interactions_Chain_B_Residues.txt "$filename"_Residues.txt > "$filename"_B_Residues_Overlapping.txt

#Save the number of unqiue residues in the reference's Chain A
Chain_A_Residue_Count="$(wc -l 2z7x_TLR1_2_Complex_Interactions_Chain_A_Residues.txt | cut -f 1 -d ' ')"

#Save the number of unqiue residues in the reference's Chain B
Chain_B_Residue_Count="$(wc -l 2z7x_TLR1_2_Complex_Interactions_Chain_B_Residues.txt | cut -f 1 -d ' ')"

#Save the number of query resiude sthat overlapped with reference's Chain A
Query_Chain_A_Count="$(wc -l "$filename"_A_Residues_Overlapping.txt | cut -f 1 -d ' ')"

#Save the number of query resiude sthat overlapped with reference's Chain B
Query_Chain_B_Count="$(wc -l "$filename"_B_Residues_Overlapping.txt | cut -f 1 -d ' ')"

#Calculate the percentage of reference's A Chain residues that overlap wtih the query's residues
Chain_A_Percentage=$(echo "scale=2;$Query_Chain_A_Count*100/$Chain_A_Residue_Count" |bc)

#Calculate the percentage of reference's B Chain residues that overlap wtih the query's residues
Chain_B_Percentage=$(echo "scale=2;$Query_Chain_B_Count*100/$Chain_B_Residue_Count" |bc)

#Save the filename and percetnage values in a file
echo $filename $Chain_A_Percentage $Chain_B_Percentage >> Docking_Percentages.txt
done

##Expanded Output File Versions
for file in 2z7x_TLR1_2_Complex_Any_Model_*_Interactions.txt
do filename=${file%_*}
grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" $file | while read p; do echo $p; done | cut -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > "$filename"_Residues.txt
Query_Chain_A_Total_Count="$(grep "A" "$filename"_Residues.txt | wc -l | cut -f 1 -d ' ')"
Query_Chain_B_Total_Count="$(grep "B" "$filename"_Residues.txt | wc -l | cut -f 1 -d ' ')"
grep -Fx -f ../2z7x_TLR1_2_Complex_Interactions_Chain_A_Residues.txt "$filename"_Residues.txt > "$filename"_A_Residues_Overlapping.txt
grep -Fx -f ../2z7x_TLR1_2_Complex_Interactions_Chain_B_Residues.txt "$filename"_Residues.txt > "$filename"_B_Residues_Overlapping.txt
Chain_A_Residue_Count="$(wc -l ../2z7x_TLR1_2_Complex_Interactions_Chain_A_Residues.txt | cut -f 1 -d ' ')"
Chain_B_Residue_Count="$(wc -l ../2z7x_TLR1_2_Complex_Interactions_Chain_B_Residues.txt | cut -f 1 -d ' ')"
Query_Chain_A_Overlapping_Count="$(wc -l "$filename"_A_Residues_Overlapping.txt | cut -f 1 -d ' ')"
Query_Chain_B_Overlapping_Count="$(wc -l "$filename"_B_Residues_Overlapping.txt | cut -f 1 -d ' ')"
Chain_A_Percentage=$(echo "scale=2;$Query_Chain_A_Overlapping_Count*100/$Chain_A_Residue_Count" |bc)
Chain_B_Percentage=$(echo "scale=2;$Query_Chain_B_Overlapping_Count*100/$Chain_B_Residue_Count" |bc)
echo $filename $Query_Chain_A_Total_Count $Query_Chain_A_Overlapping_Count $Chain_A_Residue_Count $Chain_A_Percentage $Query_Chain_B_Total_Count $Query_Chain_B_Overlapping_Count  $Chain_B_Residue_Count $Chain_B_Percentage >> 2z7x_TLR1_2_Complex_Any_Loc_Docking_Percentages_V2.txt
done