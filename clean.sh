#!/bin/bash
source ~/.bash_profile
# FOR all file in todo, done, log older than 30min rm
for f in `find ${SHINYSNP_DIR}/logs/ ${SHINYSNP_DIR}/todo/ ${SHINYSNP_DIR}/done/ -name "*_*" -mmin +30`
do 
	rm $f
done
