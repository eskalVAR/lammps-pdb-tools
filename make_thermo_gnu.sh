#!/bin/bash
sed -n $(grep -n TotEng log.lammps | cut -f1 -d:),999999999p log.lammps > thermo_data

