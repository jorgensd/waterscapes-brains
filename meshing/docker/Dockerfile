
FROM condaforge/mambaforge 

ENV DEBIAN_FRONTEND=noninteractive

# Upgrade mamba
RUN mamba upgrade -y mamba

# Copy environment and requirements files into docker env
COPY environment.yml .

# Update environment file with new environment name
RUN mamba env update --file environment.yml --name waterscales-brains

# Install freesurfer
RUN wget https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.2.0/freesurfer_7.2.0_amd64.deb
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y language-pack-en binutils libx11-dev gettext xterm x11-apps make csh tcsh file bc xorg xorg-dev libncurses5 libjpeg62
RUN dpkg -i freesurfer_7.2.0_amd64.deb
SHELL ["mamba", "run", "-n", "waterscales-brains", "/bin/bash", "-c"]
RUN echo "source activate waterscales-brains" > ~/.bashrc

COPY ./docker/.license /usr/local/freesurfer/7.2.0/.license
ENV PATH="/usr/local/freesurfer/7.2.0/bin:${PATH}"
ENV FREESURFER_HOME="/usr/local/freesurfer/7.2.0"