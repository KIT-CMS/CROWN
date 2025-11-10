FROM rootproject/root:6.34.00-ubuntu24.04

# Install necessary packages
RUN apt-get update -y
RUN apt-get install -y ca-certificates gpg wget git make python3-pip openmpi-bin libopenmpi-dev libboost-all-dev openssh-client nlohmann-json3-dev

RUN python3 -m pip install GitPython --break-system-packages && \
    python3 -m pip install git+https://github.com/cms-nanoAOD/correctionlib.git --break-system-packages

# update cmake
RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null && \
echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ jammy main' | tee /etc/apt/sources.list.d/kitware.list >/dev/null && \
apt-get update -y && \
apt-get install -y cmake

# Install onnxruntime
RUN git clone --recursive https://github.com/microsoft/onnxruntime.git
WORKDIR /opt/onnxruntime
RUN ./build.sh --config Release --update --build --build_shared_lib --parallel --allow_running_as_root
WORKDIR /

# Install tini
RUN apt-get install -y tini

ENTRYPOINT ["/usr/bin/tini", "--"]
