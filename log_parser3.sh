#!/bin/bash
# Lammps log file thermo data parser to CSV format
# Copyright (c) 2022, Eryk Skalinski
# eskalinski@protonmail.com

echo "--------------------------------------------------------"
printf "Did you know lammps support yaml format?\n"
printf "See : https://docs.lammps.org/thermo_style.html\n"
printf "See : https://docs.lammps.org/Howto_structured_data.html\n"
echo "--------------------------------------------------------"

removeMinimizationStages=1
verbose=0
numberOfStages=0

while getopts 'i:o:hmv' OPTION; do
	case "$OPTION" in
	h)
		echo "script usage: $(basename \$0) [-h] [-m] [-v] [-i lammpslog] [-o outpufilename]" >&2
		echo "outputs : outputfilename{stage_index}.csv" >&2
		exit 1
		;;
	m)
        	echo "Keeping minimization stages"
		removeMinimizationStages=0
		;;
	i)
		logFileName="$OPTARG"
		;;
	o)
		outputFileName="$OPTARG"
		;;
	v)
		verbose=1
		;;
	?)
		echo "script usage: $(basename \$0) [-h] [-m] [-v] [-i lammpslog] [-o outpufilename]" >&2
		echo "outputs : outputfilename{stage_index}.csv" >&2
		exit 1
		;;
	esac
done

if [[ -z $logFileName ]]; then echo "ERROR: No log file specified." >&2; exit 1; fi
if [[ -z $outputFileName ]]; then echo "ERROR: No output file specified.">&2; exit 1; fi
# Find thermo_style used
thermoStyle=$(grep -i "thermo_style" "$logFileName")

if [[ $(echo "$thermoStyle" | awk '{print $2}') != "custom" ]]; then echo "ERROR: Unsupported thermo_style." >&2 ; exit 1; fi

thermoHeaders=($(echo "$thermoStyle" | awk '{for(i=3;i<=NF-1;i++) printf $i" "; print ""}'))
numberofThermoHeaders=${#thermoHeader[@]}

if [[ $verbose -eq 1 ]]; then
	for header in "${thermoHeaders[@]}"
	do
		printf "Found : %s\n" "$header"
	done
fi
if [[ $removeMinimizationStages -eq 1 ]]; then
	if [[ $verbose -eq 1 ]]; then echo "Removing minimization stages."; fi
	cat "$logFileName" | awk '/^minimize/{flag=0}flag;/^Minimization stats|^LAMMPS/{flag=1}' > "$outputFileName"
	if [[ $verbose -eq 1 ]]; then echo "Removed minimization stages."; fi
fi

if [[ $verbose -eq 1 ]]; then echo "Identifying stages."; fi
numberOfStages=$(awk '{print $1}' < "$outputFileName" |  grep -ic "^${thermoHeaders[0]}")
if [[ $verbose -eq 1 ]]; then printf "Found %s stage(s).\n" "$numberOfStages"; fi
cat "$outputFileName" | awk -v outfile="$outputFileName" -v numberofFields="$numberofThermoHeaders" '
	BEGIN { indy=0 ; numberofFields=2; OFS="," } 
	/Step/{ flag=1; indy=indy+1; next; }
	{fn=outfile indy ".csv"}
	/Loop/{ flag=0 } 
	{if (flag == 1) 
		for (i=1; i<=numberofFields; i++)
			printf "%s%s", $i, (i<numberofFields ? OFS : ORS) > fn
			
	}'
