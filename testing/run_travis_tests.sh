## Due to Travis failing on runtime, I've had to change the db over from arg-annot to groot-core-db, consequently the test reads (derived from arg-annot) won't align as well

#!/bin/bash

# install the software
go build

# download the database
./groot get -d groot-core-db -o ./db

# index the example
./groot index -p 1 -i ./db/groot-core-db.90 -o test-index

# align the test reads
./groot align -p 1 -i test-index -f testing/full-argannot-perfect-reads-small.fq.gz > out.bam

# report
./groot report -i out.bam