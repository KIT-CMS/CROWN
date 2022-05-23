#!/bin/bash
function overridePar() {
    file=$1
    arg=$2
    val=$3
    if [[ "" = $( grep -e "$arg" $file ) ]]; then
        logerror $arg is not in $file
        return 2
    fi
    if [[ $( grep -e "(^|\s|--)$arg" $file | grep -c . ) >1 ]]; then
        logerror $arg appears multiple times in $file, aborting.
        return 2
    fi
    sed -i -E "s@(^|\s|--)$arg(\s+|=)\w+@\1$arg\2$val@" $file
}
function getPar() {
    file=$1
    arg=$2
    lines=$( grep -E "(^|\s|--)$arg" $file  )
    nlines=$( echo $lines | grep -c . )

    if [[ $nlines == 0 ]]; then
        logerror $arg is not in $file
        return 2
    fi
    if [[ $nlines >1 ]]; then
        logerror $arg appears multiple times in $file, aborting.
        return 2
    fi
    echo $lines | sed -E "s@(^|\s|--)($arg)(\s+|=)(\w+).*@\4@"
}

function updateSymlink() {
    source=$1
    target=$2
    loginfo "Updating Link: $target -> $source"
    if [[ -d $source ]]; then
        [[ -L $target ]] && rm $target
        [[ -d $target ]] && rmdir $target
        ln -sn $source $target
    elif [[ -f $source ]]; then
        [[ -f $target || -L $target ]]  && rm $target
        ln -s $source $target
    else
        logerror "Invalid source: $source"
        return 1
    fi

}
function myCpuUsage() {
    top -b -n 1 -u $USER | awk 'NR>7 { sum += $9; } END { printf "%3.0f\n", sum/100; }'
}

function condwait(){
    sleep .1
    if [[ $( myCpuUsage ) -gt $( recommendCPUs ) ]]; then
        wait
    fi
}

function recommendCPUs() {
    ## recommeds % of free cpus
    avUsage=$(grep 'cpu ' /proc/stat | awk '{usage=($2+$4)*100/($2+$4+$5)} END {print usage}')
    ncpus=$(($(grep 'cpu' /proc/stat | wc -l)-1))
    # based on percentare of currently used cpus
    #echo $avUsage $ncpus | awk '{print int((1-$1/100)*$2*.5)}'
    # fixed to 25%
    echo $ncpus | awk '{print int($1/4)}'
}


function capargs() {
    outstr=""
    for arg in $( echo $@ | cut -d' ' -f1-4 ); do
        if ispath $arg; then
            arg=$(shortenpath $arg )
        fi
        if [[ ${#arg} -gt 30 ]]; then
            outstr=$outstr+${arg:0:30}
        else
            outstr=$outstr+$arg
        fi
    done
    echo $outstr | sed 's@\./@@' | sed -E "s@/@:@g" | sed -E "s@^[-:+]+@@g"  #| sed -E "s@ +@+@" 
}
function ispath() {
    if [[ -d $@ || -f $@ ]]; then
        return 0
    else
        return 1
    fi
}
function shortenpath() (
    [[ -n $ZSH_VERSION ]] && setopt +o nomatch
    pathstr=""
    if [[ $1 = /* ]]; then
        pathstr="/"
        findbase="/"
    else
        pathstr=""
        findbase="."
    fi
    for level in $( echo $1 | sed "s@/@ @g" ); do
        pathstr=${pathstr}${level:0:3}
        level=${level:3}
        while [[ $(echo ${pathstr}* 2> /dev/null | wc -l) -gt 1 ]]; do
            pathstr=${pathstr}${level:0:1}
            level=${level:1}
        done
        [[ $level == "" ]] || pathstr=${pathstr}
        pathstr=${pathstr}/
    done
    echo $pathstr
)

# colored output for viewing logs
LESS=-R
function loginfo {
    echo -e "\e[46m[INFO]\e[0m" $( date +"%y-%m-%d %R" ): $@ | tee -a $( pwd )/output/log/logandrun/event.log
}
function logwarn {
    echo -e "\e[45m[WARN]\e[0m" $( date +"%y-%m-%d %R" ): $@ | tee -a $( pwd )/output/log/logandrun/event.log
}
function logerrormsg {
    echo -e "\e[41m[ERROR]\e[0m" $( date +"%y-%m-%d %R" ): $@ | tee -a $( pwd )/output/log/logandrun/event.log
}
function logerror {
    logerrormsg $@
    return 1
}
function logandrun() (
    # set bash to exit if as soon as a part of a pipe has a non-zero exit value
    set -o pipefail
    set +e
    # check if this is the top level logandrun script, so multiple notifications can be supressed
    if [[ -z $LOGANDRUN_TOPLEVELSET ]]; then
        export LOGANDRUN_TOPLEVELSET=1
        LOGANDRUN_IS_TOP_LEVEL=1
    fi
    # set the name of the logfile based on the command
    logfile=$( pwd )/output/log/logandrun/$( capargs "$@" ).log
    # print the startmessage and log it
    echo -e "\e[43m[RUN]\e[0m" $( date +"%y-%m-%d %R" ): $@ | tee -a $( pwd )/output/log/logandrun/event.log | tee -a $logfile
    # evaluate the current date in seconds as the start date
    start=`date +%s`
    #######
    # execute the command and log it
    (    
    set -e
    $@ 2>&1 |  tee -a $logfile
    )
    # capture the return code ( without  pipefail this would be the exit code of tee )
    return_code=$?
    end=`date +%s`
    # evaluate the end data
    # if the there was no error...
    if [[ $return_code == 0 ]]; then
        # print the message and log it
        echo -e "\e[42m[COMPLETE]\e[0m" $( date +"%y-%m-%d %R" ): $@ "     \e[104m{$((end-start))s}\e[0m" | tee -a $( pwd )/output/log/logandrun/event.log | tee -a $logfile
    else
        # print a message with the return code
        logerrormsg Error Code $return_code  $@ "     \e[104m{$((end-start))s}\e[0m"  | tee -a $( pwd )/output/log/logandrun/event.log | tee -a $logfile
    fi
    # check if the script has been running for more than 2400s = 40 min AND
    # there is a script to notify the user
    if [[ $((end-start)) -gt 2400 && $LOGANDRUN_IS_TOP_LEVEL==1 ]] && hash sendmsg.py 2>/dev/null ; then
        #notify the user
        sendmsg.py "$(hostname) $( date +"%y-%m-%d %R" ) {$((end-start))s} $return_code: $@ " 2>/dev/null
    fi
    if [[ $LOGANDRUN_IS_TOP_LEVEL==1 ]]; then
        unset LOGANDRUN_TOPLEVELSET
        unset LOGANDRUN_IS_TOP_LEVEL
    fi
    return $return_code
)
function logclean () {
    find output/log -type f -iname "*.log" -delete
}

# this function makes sure all the output directories exit
function ensureoutdirs() {
    [[ -d output ]] || mkdir output
    pushd output  >/dev/null
    for folder in datacards log log/logandrun plots  shapes ml signalStrength; do
        [[ ! -d $folder ]] && mkdir $folder
    done
    [[ -d log/condorShapes ]] || mkdir log/condorShapes
    popd  >/dev/null
}
