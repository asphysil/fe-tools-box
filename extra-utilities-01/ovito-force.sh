#!/bin/bash

#rm dump.xyz 

#cp bfo_dmi-0.1.xyz dump.xyz 

sed -i '' -e 's/c_outsp\[1\]/fx/g' dump.xyz 
sed -i '' -e 's/c_outsp\[2\]/fy/g' dump.xyz
sed -i '' -e 's/c_outsp\[3\]/fz/g' dump.xyz 
