# start from rootproject/root:6.30.04-ubuntu22.04 to be changed to rootproject/root:6.32.02-ubuntu22.04
FROM rootproject/root:6.36.00-ubuntu25.04

# Install necessary packages
RUN apt-get update -y
RUN apt-get install -y ca-certificates gpg wget git make python3-pip openmpi-bin libopenmpi-dev libboost-all-dev openssh-client nlohmann-json3-dev

RUN python3 -m pip install GitPython && python3 -m pip install git+https://github.com/cms-nanoAOD/correctionlib.git

# update cmake
RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null && \
echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ jammy main' | tee /etc/apt/sources.list.d/kitware.list >/dev/null && \
apt-get update -y && \
apt-get install -y cmake

# Install onnxruntime
RUN git clone --recursive https://github.com/microsoft/onnxruntime.git
RUN cd onnxruntime && ./build.sh --config Release --update --build --build_shared_lib --parallel --allow_running_as_root

# Install tini
RUN apt-get install -y tini

ENTRYPOINT ["/usr/bin/tini", "--"]
