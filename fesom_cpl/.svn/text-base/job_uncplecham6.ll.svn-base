#!/bin/bash
# @ job_name         = fesom_uncpl
# @ output           = $(job_name).out
# @ error            = $(job_name).err
# @ environment      = COPY_ALL
# @ job_type         = parallel
# @ notification     = never ##always 
# @ node             = 2
# @ tasks_per_node   = 64
# @ resources        = ConsumableMemory(300Mb)
# @ network.mpi      = sn_all,not_shared,us
# @ task_affinity    = cpu(1)
# @ class            = express
# @ wall_clock_limit = 00:20:00
# @ account_no       = ab0046
# @ queue

# ======================================= USER SECTION ========================================


export VT_MODE=STAT:TRACE
export VT_FILE_PREFIX=fesom_uncpl
export VT_VERBOSE=2
export VT_MAX_FLUSHES=0
export VT_BUFFER_SIZE=2G


export runid=fesom	#FS200 hier ist doch vieles Weiteres auch unnoetig, oder???

echo $runid

export codepath=${HOME}/fesom_echam6_oasis3-mct/fesom_cpl
export resultpath=/work/${GROUP}/${USER}/fesom_uncplecham6/
echo 'welcome on'
hostname

echo $codepath

echo 'do a full year'
cd $codepath

#cp inputfile_10d inputfile
rm -f goodfile

poe ./fesom.x			#invokes Parallel Operating Environment (POE)

export ex=`find . -name goodfile -print`

if [ ${#ex} -eq 0 ];  then
 echo 'Oh wie schade. Hier hat etwas nicht funktioniert.'
 exit
fi

echo 'Start next job'



exit
