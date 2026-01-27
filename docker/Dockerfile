FROM continuumio/miniconda3
SHELL ["/bin/bash", "-c"]
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        cmake \
        make \
        g++ \
    && rm -rf /var/lib/apt/lists/*
RUN mkdir /tools
WORKDIR /tools
COPY environment.yaml .
RUN . /opt/conda/bin/activate && \
    conda env create -n pantax -f environment.yaml && \
    conda install rust=1.82 clangdev=19 glpk  && \
    conda clean -ya

ENV GUROBI_HOME=/opt/conda/envs/pantax
RUN git clone https://github.com/LuoGroup2023/PanTax.git -b all-rust && \
    bash install.sh && \
    cp scripts/pantax-rg scripts/pantax-utils scripts/pantax /opt/conda/envs/pantax/bin/ && \
    cp scripts/*py /opt/conda/envs/pantax/bin/
ENV PATH=/opt/conda/envs/pantax/bin/:$PATH
ENV LD_LIBRARY_PATH=/opt/conda/envs/pantax/lib:$LD_LIBRARY_PATH
CMD ["conda", "run", "-n", "pantax"]


# docker build -t pantax:v1 .
# # 1. run directly in your path with data
# docker run -v $(dirname $PWD):/mnt -w /mnt/$(basename $PWD) pantax:v1 pantax -f ../example_genomes_info.txt -s -p -r short_reads.fq.gz --species --strain -t 64 --solver highs
# # 2. start an interactive docker container session and run in your path with data
# docker run -it --rm -v $(dirname $PWD):/mnt -w /mnt/$(basename $PWD) -v /var/run/docker.sock:/var/run/docker.sock pantax:v1 /bin/bash
# conda activate pantax
# pantax -f ../example_genomes_info.txt -s -p -r short_reads.fq.gz --species --strain -t 64 --solver highs

# gurobi does not work, it need WSL license for docker
# docker run \
#   --rm \
#   --net=host \
#   --uts=host \
#   --ipc=host \
#   --pid=host \
#   -v $(dirname $PWD):/mnt \
#   -v /home/wenhai/application/gurobi1100/linux64/license_file/002:/opt/gurobi \
#   -e GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic \
#   -e GRB_MODEL=ENV \
#   -w /mnt/$(basename $PWD) \
#   pantax:v1 \
#   pantax -f ../example_genomes_info.txt -s -p -r short_reads.fq.gz --species --strain -t 64
