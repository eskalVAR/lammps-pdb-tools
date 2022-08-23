#!/bin/bash
# --------------------------
# Eryk Skalinski 2022
# eskalinski@protonmail.com
# --------------------------
# VERSION 0.5
# --------------------------
# CHANGELOG
# --------------------------
# 0.5 Overhauled internal report of data file
# 0.4 Added counter-ion generation and addition to data file
# 0.3 Added flags and overall rebasing
# 0.2 Improved missing co-efficient discovery
# 0.1 Intial version
# --------------------------
# Dependencies:
# Perl
# charmm2lammps.pl (Available from lammps/tools/ch2lmp)
# --------------------------
# Other requirements:
# Force-field files
# - par_*.prm
# - top_*.rtf
# *.pdb, *.psf, and *.crd files generated from CHARMM-GUI PDB Reader
# (Available at www.charmm-gui.org)
needReRun=0
fix_lammps_data(){
	sed -i "2 i\ " "$1"
	sed -i "8 i\ " "$1"
	sed -i "14 i\ " "$1"
	sed -i -E -e '/Masses|Atoms|Bond Coeffs|Bonds|Angle Coeffs|Angles|Dihedral Coeffs|Dihedrals|Improper Coeffs|Impropers/G' -e '/Masses|Atoms|Bond Coeffs|Bonds|Angle Coeffs|Angles|Dihedral Coeffs|Dihedrals|Improper Coeffs|Impropers/i\ ' "$1"
}
insert_lammps_mass(){
	coeff_start=$(grep -n "$2" "$1"| cut -f1 -d:)
	coeff_start=$((coeff_start))
	coeff_end=$(grep -n "$3" "$1" | cut -f1 -d:)
	coeff_end=$((coeff_end-1))
	final_coeff=$(sed -n "$coeff_end"p "$1" | awk '{print $1}')
	final_coeff=$(($final_coeff + 1))
	echo $final_coeff
	cropped_line=$(echo "$4" | awk '{printf "%10s%3s %s",$2,$3,$4}')
	new_line=$(printf "%7s %s" "$final_coeff" "$cropped_line")
	coeff_end=$((coeff_end+1))
	sed -i "$coeff_end i\ $new_line" "$1"
}
insert_lammps_pair_coeff(){
	coeff_start=$(grep -n "$2" "$1"| cut -f1 -d:)
	coeff_start=$((coeff_start))
	coeff_end=$(grep -n "$3" "$1" | cut -f1 -d:)
	coeff_end=$((coeff_end-1))
	final_coeff=$(sed -n "$coeff_end"p ./*.data | awk '{print $1}')
	final_coeff=$(($final_coeff + 1))
	echo $final_coeff
	cropped_line=$(echo "$4" | awk '{printf "%16s%17s%17s%17s%3s %s",$2,$3,$4,$5,$6,$7}')
	new_line=$(printf "%7s %s" "$final_coeff" "$cropped_line")
	coeff_end=$((coeff_end+1))
	sed -i "$coeff_end i\ $new_line" "$1"
}

insert_lammps_atom(){

	coeff_start=$(grep -n "$2" "$1"| cut -f1 -d:)
	coeff_start=$((coeff_start))
	coeff_end=$(grep -n "$3" "$1" | cut -f1 -d:)
	coeff_end=$((coeff_end-1))
	atom_coeff=$(sed -n "$coeff_end"p ./*.data | awk '{print $1}')
	atom_coeff=$(($atom_coeff + 1))
	molecule_coeff=$(sed -n "$coeff_end"p ./*.data | awk '{print $2}')
	molecule_coeff=$(($molecule_coeff + 1))
	cropped_line=$(echo "$4" | awk '{printf "%10s%17s%17s%17s%3s %s",$4,$5,$6,$7,$8,$9}')
	new_line=$(printf "%7s%8s%6s%s" "$atom_coeff" "$molecule_coeff" "$5" "$cropped_line")
	coeff_end=$((coeff_end+1))
	sed -i "$coeff_end i\ $new_line" "$1"
}
pre_proc(){
	link_ch2lmp_files(){
		ln -f ../charmm2lammps.pl . || exit
		ln -f ../*$forceFieldName* . || exit
		ln -f "../${projectName}_pdbreader.crd" . || exit
		ln -f "../${projectName}_pdbreader.psf" . || exit
	}
	clean_link_files(){
		rm charmm2lammps.pl 
		rm ./*$forceFieldName*
		rm "./${projectName}_pdbreader.crd"
		rm "./${projectName}_pdbreader.psf"
	}
	# Uncompressing
	printf 'Uncompressing charmm-gui.tgz\n'
	tar zxf ./charmm-gui.tgz
	printf 'Renaming folder contents\n'
	cd ./charmm-gui-*/. || exit 1
	cp step1_pdbreader.pdb "../${projectName}_pdbreader.pdb"
	cp step1_pdbreader.psf "../${projectName}_pdbreader.psf"
	cp step1_pdbreader.crd "../${projectName}_pdbreader.crd"
	cd ..

	# No Ions
	mkdir no_ions
	cd no_ions || exit 1
	printf 'Hard sym-linking files\n'
	link_ch2lmp_files "$projectName"
	printf 'Attempting to convert %s to lammps data\n' "$projectName"
	printf '########################################\n'
	c2l_out=$(perl ./charmm2lammps.pl "$forceFieldName" "${projectName}_pdbreader")
	printf 'charmm2lammps.pl reports the following:\n'
	printf '########################################\n'
	echo "$c2l_out" | grep "was not found"
	printf '########################################\n'
	printf 'Cleaning sym-linked files\n'	
	printf '########################################\n'
	clean_link_files "$projectName"
	cd ..
	cp ./no_ions/*ctrl.psf "./${projectName}_ctrl.psf"
	grep -v '^$' ./no_ions/*.data > "./${projectName}_ctrl.data"
	cp "./${projectName}_ctrl.data" "./${projectName}.data"
	# Ions
	if [[ $ions -eq 1 ]]; then
		printf 'Generating counter ions\n'
		printf '########################################\n'
		mkdir water
		cd water || exit
		link_ch2lmp_files "$projectName"
		perl ./charmm2lammps.pl "$forceFieldName" "${projectName}_pdbreader" -water -ions > /dev/null 2> /dev/null
		clean_link_files "$projectName"
		cd ..
		mkdir ions
		cd ions || exit
		cp "../${projectName}_ctrl.data" ions.data
		
		OLDIFS=$IFS
		IFS="
		"

		NaIons=($(grep "SOD" ../water/*.data))
		ClIons=($(grep "CLA" ../water/*.data))
		NaTotal=$((${#NaIons[@]}-2))
		ClTotal=$((${#ClIons[@]}-2))
		NaAtoms=("${NaIons[@]:2}")
		ClAtoms=("${ClIons[@]:2}")
		IFS=$OLDIFS
		printf 'Adding %d Na ions and %d Cl Ions\n' "$NaTotal" "$ClTotal"
		NaID=$(insert_lammps_mass "./ions.data" "Masses" "Pair Coeffs" "${NaIons[0]}")
		ClID=$(insert_lammps_mass "./ions.data" "Masses" "Pair Coeffs" "${ClIons[0]}")
		junk=$(insert_lammps_pair_coeff "./ions.data" "Pair Coeffs" "Atoms" "${NaIons[1]}")
		junk=$(insert_lammps_pair_coeff "./ions.data" "Pair Coeffs" "Atoms" "${ClIons[1]}")
		if [[ ${#NaAtoms[@]} -ne 0 ]]; then
			for atom in "${NaAtoms[@]}"
			do
				insert_lammps_atom "./ions.data" "Atoms" "Bond Coeffs" "$atom" "$NaID"
			done
		fi
		if [[ ${#ClAtoms[@]} -ne 0 ]]; then
			for atom in "${ClAtoms[@]}"
			do
				insert_lammps_atom "./ions.data" "Atoms" "Bond Coeffs" "$atom" "$ClID"
			done
		fi

		printf 'Insertion finished and header ammended\n'
		printf '########################################\n'
		mv ./ions.data "../${projectName}.data"
		cd ..
	fi
}
internal_report(){
	coeff_start=$(grep -n "$2" "$1" | cut -f1 -d:)
	coeff_start=$((coeff_start+1))
	coeff_end=$( grep -n "$3" "$1" | cut -f1 -d:)
	coeff_end=$((coeff_end-1))
	value_start=$((coeff_end+2))
	value_end=$(grep -n "$4" "$1" | cut -f1 -d:)
	value_end=$((value_end-1))
	printf "%s:\n" "$2"
	coeffs=$(sed -n "$coeff_start","$coeff_end"p "$1" | awk '{print $1}')
	values=$(sed -n "$value_start","$value_end"p "$1" | awk '{print $2}')
	required=$(comm -23 <(sort -u <<<"${values[*]}") <(sort -u <<<"${coeffs[*]}"))
	dummy=$(echo "${coeffs[@]}" | awk '{for(i=p+1; i<$1; i++) print i} {p=$1}' )
	for i in ${dummy[@]}
	do
	match=0
		for j in ${required[@]}
		do
		if [ "${i}" == "${j}" ]
		then
			match=1
			break
		fi
		done
	if [ "${match}" == 0 ]
	then
		filt_dummy=(${filt_dummy[@]} $i)
	fi
	done
	for req in ${required[@]}
	do
		fixed_req=(${fixed_req[@]} $req)
	done
	if [[ -n "${filt_dummy}" ]]; then
		needReRun=1
		for value in "${filt_dummy[@]}"
		do
			printf '%d can be made into a dummy\n' "$value"
		done
	fi
	if [[ -n "${fixed_req[@]}" ]]; then
		needReRun=1
		for value in "${fixed_req[@]}"
		do
			printf '%d requires a declaration\n' "$value"
		done
	fi
}
total_report(){
	printf 'Internal check reports the following:\n'
	printf '(Run with lammps/ovito to double check)\n'
	printf '########################################\n'
	internal_report "$projectName.data" "Bond Coeffs" "Bonds" "Angle Coeffs"
	printf '########################################\n'
	internal_report "$projectName.data" "Angle Coeffs" "Angles" "Dihedral Coeffs"
	printf '########################################\n'
	internal_report "$projectName.data" "Dihedral Coeffs" "Dihedrals" "Improper Coeffs"
	printf '########################################\n'
	if [[ $needReRun -eq 1 ]]; then
		printf 'Fix Coefficients before continuing\n'
		printf 'Exiting\n'
		exit 1
	else
		printf 'All Coefficients appear valid\n'
		printf 'Continuing\n'
	fi
}
finalMassID=0
extract_pair_coeffs(){
	start=$(grep -n "Pair Coeffs" "$1" | cut -f1 -d:)
	end=$(grep -n "Atoms" "$1" | cut -f1 -d:)
	end=$((end - 1))
	# Extract
	# Delete from file
	start=$((start + 1))
	if [[ $goModel -eq 1 ]]; then
		sed -n "$start","$end"p "$1" | awk -v HID="$finalMassID" > "${projectName}_parameters_go.txt" '
		BEGIN { NID=HID+1; OID=HID+2 }
		{printf "pair_coeff%8s%8slj/charmm/col/charmm/implicit%8s%17s%17s%17s%3s %s%s",$1,$1,$2,$3,$4,$5,$6,$7,ORS}
		END {
		printf "##### hbond/dreiding/lg for GO-Model Hydrogen Bonds%s",ORS
		printf "%10s%8s%8s%23s%8s%s", "#", "O", "N", "#", "H", ORS
		printf "pair_coeff%8s%8s      hbond/dreiding/lj%8s%8s%8s%8s%8s%8s%8s%8s%s",OID,NID,HID, "j", "${eps}", "3.1626906", "4", "9.0","11.0", "90.0", ORS
		}
		'
	else
		sed -n "$start","$end"p "$1" | awk > "${projectName}_parameters.txt" '
		{printf"pair_coeff%8s%8s%8s%17s%17s%17s%3s %s%s",$1,$1,$2,$3,$4,$5,$6,$7,ORS}'
	fi
	start=$((start - 1))
	sed -i -e "${start},${end}d" "$1"
}
handle_hbonds(){
	# For each reisdue:
	# First atom is the N
	# Second atom is the H
	# Last atom is the O
	# H,N,O
	finalMass=$(grep -n "Atoms" "$1" | cut -f1 -d:)
	finalMass=$((finalMass - 1))

	finalMassID=$(sed -n "$finalMass"p "$1" | awk '{print $1}')
	finalMassID=$((finalMassID + 1))
	tmpfile=$(mktemp /tmp/charmm2lammpsgomodel.XXXXXX)

	start=$(grep -n "Atoms" "$1" | cut -f1 -d:)
	end=$(grep -n "Bond Coeffs" "$1" | cut -f1 -d:)
	start=$((start+1))
	end=$((end-1))
	# Substitue IDs
	sed -n "$start","$end"p "$1" | awk -v HID="$finalMassID" -v aaS="$aminoAcids" > "$tmpfile" '
	BEGIN { chainID=0; aminoAcidID=0; fixnext=0; NID=HID+1; OID=HID+2}
	{
		if ($9 != "SOD" && $9 != "CLA")
		{
			if (fixnext == 1)
			{
				fixnext=0
				printf "%8s%8s%6s", $1, $2, HID
				printf "%10s%17s%17s%17s%3s %s%s",$4,$5,$6,$7,$8,$9,ORS
			}
			else if ($2 == aa)
			{
				if($9 == "O")
				{
					printf "%8s%8s%6s", $1, $2, OID
					printf "%10s%17s%17s%17s%3s %s%s",$4,$5,$6,$7,$8,$9,ORS
					# NID=NID+3
					# HID=HID+3
					# OID=OID+3
				}
				else
				{
					printf "%8s%8s%6s%10s%17s%17s%17s%3s %s%s",$1,$2,$3,$4,$5,$6,$7,$8,$9,ORS
				}
			}
			else
			{
				printf "%8s%8s%6s", $1, $2, NID
				printf "%10s%17s%17s%17s%3s %s%s",$4,$5,$6,$7,$8,$9,ORS
				aa=$2
				fixnext=1
			}
		}
		else
		{
			printf "%8s%8s%6s%10s%17s%17s%17s%3s %s%s",$1,$2,$3,$4,$5,$6,$7,$8,$9,ORS
		}
	}
	'
	# Merge Files
	section_start=$(grep -n "Atoms" "$1" | cut -f1 -d:)
	section_start=$((section_start + 1))
	section_end=$( grep -n "Bond Coeffs" "$1" | cut -f1 -d:)
	section_end=$((section_end - 1))
	sed -i "$section_start","$section_end"d "$1"
	section_start=$((section_start - 1))
	sed -i -e "${section_start}r $tmpfile" "$1"

	# Add Masses
	junk=$(insert_lammps_mass "$1" "Masses" "Pair Coeffs" "0 1.008 # H")
	junk=$(insert_lammps_mass "$1" "Masses" "Pair Coeffs" "0 14.007 # N")
	junk=$(insert_lammps_mass "$1" "Masses" "Pair Coeffs" "0 15.999 # O")

	rm "$tmpfile"
}
usage() {
	echo "script usage: $(basename \$0) [-c] [-h] [-g] [-i] [-p {project name}] [-f {forcefield name}] [-r]" >&2
	exit 1
}
reportOnly=0
ions=0
goModel=0
while getopts 'p:f:gihrc' OPTION; do
	case "$OPTION" in
	h)
		# Help
		usage
		;;
	c)
		# Clean
		printf "WARNING CLEANING\n"
		printf "WHAT WILL BE REMOVED:\n"
		printf "ALL SUB DIRECTORIES\n"
		printf "*.crd\n"
		printf "*.pdb\n"
		printf "*.psf\n"
		printf "*.data\n"
		printf "*.txt\n"
		read -p "Are you sure? " -n 1 -r
		echo 
		if [[ $REPLY =~ ^[Yy]$ ]]
		then
			rm ./*.crd 2> /dev/null
			rm ./*.pdb 2> /dev/null
			rm ./*.psf 2> /dev/null
			rm ./*.data 2> /dev/null
			rm ./*.txt 2> /dev/null
			rm -rf ./*/ 2> /dev/null
		else
			exit 0
		fi
		;;
	g)
		# Basic equal probability Go model
		goModel=1
		;;
	i)
		# Ions
		ions=1
		;;
	p)
		# Project name
		projectName="$OPTARG"
		;;
	f)
		# Force field
		forceFieldName="$OPTARG"
		;;
	r)
		# Just report
		reportOnly=1
		printf '-r supplied, only running report\n'
		;;
	?)
		usage
		;;
	esac
done
set -eo pipefail
if [ -z "$forceFieldName" ]
then
	printf 'No force field supplied with -f flag\n'
	printf 'Exiting\n'
	exit 1
	
fi

if [ -z "$projectName" ]
then
	printf 'No project name supplied with -p flag\n'
	printf 'Exiting\n'
	exit 1
	
fi
if [[ $reportOnly -eq 1 ]]; then
	total_report
	exit 0
fi

if [ -f "${projectName}_ctrl.data" ] && [ -f "${projectName}.data" ] && [ -f "${projectName}_ctrl.psf" ]; then
	printf 'Found existing files, skipping pre-processing\n'
else
	pre_proc "$projectName"
fi

total_report


if [[ $goModel -eq 1 ]]; then
	printf 'Applying Go model\n'
	handle_hbonds "$projectName.data"
fi

printf 'Generating pair coefficients file\n'
extract_pair_coeffs "$projectName.data"

printf 'Reformatting data file\n'
fix_lammps_data "$projectName.data"

printf 'Done\n'

