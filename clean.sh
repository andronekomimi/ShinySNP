#!/bin/bash

# FOR all file in todo, done, log older than 30min rm

for f in `find logs/ todo/ done/ -name "*_*" -mmin +30`
do 
	rm $f
done
