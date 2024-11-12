
FROM condaforge/mambaforge:24.9.0-0

ENV DEBIAN_FRONTEND=noninteractive
ARG FREESURFER_VERSION=7.4.0


# Upgrade mamba
RUN mamba upgrade -y mamba

# Copy environment and requirements files into docker env
COPY environment.yml .
# Update environment file with new environment name
RUN mamba env update --file environment.yml --name waterscales-brains

RUN apt-get update && apt-get install -y language-pack-en binutils libx11-dev gettext xterm x11-apps make csh tcsh file bc xorg xorg-dev libncurses5 libjpeg62 lsb-release 
# Install freesurfer
RUN wget https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/${FREESURFER_VERSION}/freesurfer_ubuntu$(lsb_release -sr | cut -d "." -f 1)-${FREESURFER_VERSION}_amd64.deb

# Install freesurfer specific packages
RUN apt-get install -y libquadmath0 libwayland-cursor0 libxcb-icccm4 libxcb-image0 libxcb-render-util0 libxcb-render0 libxcb-xinerama0  libxcb-xinput0 libxcb-xkb1 libxkbcommon-x11-0 libxkbcommon0 libxcb-keysyms1

ENV DEBIAN_FRONTEND=noninteractive
RUN dpkg -i freesurfer_ubuntu$(lsb_release -sr | cut -d "." -f 1)-${FREESURFER_VERSION}_amd64.deb
SHELL ["mamba", "run", "-n", "waterscales-brains", "/bin/bash", "-c"]
RUN echo "source activate waterscales-brains" > ~/.bashrc

COPY ./docker/.license /usr/local/freesurfer/${FREESURFER_VERSION}/.license
ENV PATH="/usr/local/freesurfer/${FREESURFER_VERSION}/bin:${PATH}"
ENV FREESURFER_HOME="/usr/local/freesurfer/${FREESURFER_VERSION}"

RUN apt -y install mesa-utils
ENV MESA_LOADER_DRIVER_OVERRIDE=i965
