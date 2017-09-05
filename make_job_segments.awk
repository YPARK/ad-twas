#!/usr/bin/awk -f
BEGIN {
    curr = 0
    for(j = 0; j < (NTOT - CHUNK); j += CHUNK) {
	curr = j + CHUNK
	print (j + 1) ":" curr
    }
    if(curr < NTOT)
	print (curr + 1) ":" NTOT
}
