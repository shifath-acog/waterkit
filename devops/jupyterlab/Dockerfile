FROM nvidia/cuda:11.8.0-devel-ubuntu22.04

WORKDIR /home/waterkit

RUN apt-get update && \
    apt-get install -y curl wget git python3 python3-pip build-essential autoconf automake libtool vim && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda && \
    rm -rf /var/lib/apt/lists/* && rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/miniconda/bin:$PATH"

COPY environment.yaml /home/
RUN conda tos accept --override-channels --channel conda-forge \
    --channel https://repo.anaconda.com/pkgs/main \
    --channel https://repo.anaconda.com/pkgs/r && \ 
    conda env create -f /home/environment.yaml && \
    conda clean --all --yes

RUN echo "source /miniconda/etc/profile.d/conda.sh && conda activate waterkit" >> /root/.bashrc

SHELL ["/bin/bash", "--login", "-c"]

RUN git clone https://github.com/forlilab/waterkit /home/waterkit && \
    /bin/bash --login -c "source /miniconda/etc/profile.d/conda.sh && conda activate waterkit && pip install -e /home/waterkit"

RUN cp /home/waterkit/scripts/*.py /miniconda/envs/waterkit/bin/ && \
    chmod +x /miniconda/envs/waterkit/bin/*.py

COPY wk_create_grid_protein_file.py /miniconda/envs/waterkit/bin/wk_create_grid_protein_file.py
COPY wk_prepare_receptor.py /miniconda/envs/waterkit/bin/wk_prepare_receptor.py

RUN chmod +x /miniconda/envs/waterkit/bin/wk_create_grid_protein_file.py
RUN chmod +x /miniconda/envs/waterkit/bin/wk_prepare_receptor.py

RUN mkdir -p /home/waterkit/WATERKIT

WORKDIR /home/waterkit/WATERKIT

RUN conda install -c bioconda autodock && \
    conda clean --all --yes

RUN cd /home/waterkit/autodocksuite-4.2.6-src/autogrid && \
    autoreconf -i && \
    mkdir x86_64Linux2 && \
    cd x86_64Linux2 && \
    ../configure && \
    make && \
    make install

RUN apt-get update && \
    apt-get install -y gnupg2 curl && \
    curl -s -L https://nvidia.github.io/nvidia-container-runtime/gpgkey | apt-key add - && \
    distribution=$(. /etc/os-release; echo $ID$VERSION_ID) && \
    curl -s -L https://nvidia.github.io/nvidia-container-runtime/$distribution/nvidia-container-runtime.list | tee /etc/apt/sources.list.d/nvidia-container-runtime.list && \
    apt-get update && \
    apt-get install -y nvidia-container-runtime && \
    mkdir -p /etc/OpenCL/vendors && echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd


COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

CMD ["bash", "--login", "-c", "tail -f /dev/null"]