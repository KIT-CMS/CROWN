# add ~/.local/bin to path if it is not already there
pathappend() {
  for ARG in "$@"
  do
    if [ -d "$ARG" ] && [[ ":$PATH:" != *":$ARG:"* ]]; then
        PATH="${PATH:+"$PATH:"}$ARG"
    fi
  done
}
pathappend "~/.local/bin"
if uname -a | grep -E 'el7' -q
then
    # source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-clang10-opt/setup.sh
    source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc9-opt/setup.sh
    # source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc9-dbg/setup.sh
else
    echo "You are not running on CentOS7, things will propably break..."
fi
