#!/bin/bash

source params.sh

scp ${HOST}:${HOSTBASEDIR}/constraint_strsets/*.tab ${DATADIR}/constraint_strsets/
scp ${HOST}:${HOSTBASEDIR}/constraint_strsets/strsets/* ${DATADIR}/constraint_strsets/

