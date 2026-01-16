# Dockerfile para contenedor estilo MFEM developer-cpu
# Basado en Ubuntu con herramientas de desarrollo y MFEM

FROM ubuntu:22.04

# Evitar prompts interactivos
ENV DEBIAN_FRONTEND=noninteractive

# Crear usuario no root
ARG USER=euler
ARG UID=1000
ARG GID=1000

# Actualizar e instalar dependencias básicas
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    vim \
    gfortran \
    libopenmpi-dev \
    openmpi-bin \
    libmetis-dev \
    libblas-dev \
    liblapack-dev \
    zlib1g-dev \
    libarchive-tools \
    supervisor \
    python3-pip \
    git-lfs \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Crear usuario euler
RUN groupadd -g ${GID} ${USER} && \
    useradd -m -u ${UID} -g ${GID} -s /bin/bash ${USER}

# Cambiar a root para instalación de VSCode
USER root
WORKDIR /opt

# Instalar OpenVSCode Server (opcional, para desarrollo web)
RUN mkdir -p /opt/archives
WORKDIR /opt/archives

# OpenVSCode server
RUN curl -L https://github.com/gitpod-io/openvscode-server/releases/download/openvscode-server-v1.69.1/openvscode-server-v1.69.1-linux-x64.tar.gz >  /opt/archives/openvscode-server-v1.69.1-linux-x64.tar.gz && \
  tar xzf openvscode-server-v1.69.1-linux-x64.tar.gz && chown -R euler:euler openvscode-server-v1.69.1-linux-x64

USER euler
WORKDIR /home/euler


# Extensiones de OpenVSCode
RUN mkdir -p ${HOME}/.openvscode-server/extensions
RUN cd ${HOME}/.openvscode-server/extensions && \
     curl -4 -L https://open-vsx.org/api/ms-python/python/2022.16.1/file/ms-python.python-2022.16.1.vsix > python.vsix && \
     bsdtar -xvf python.vsix extension && \
     mv extension ms-python.python-2022.16.1

RUN cd ${HOME}/.openvscode-server/extensions && \
    curl -L https://open-vsx.org/api/ms-toolsai/jupyter/2022.7.1001890528/file/ms-toolsai.jupyter-2022.7.1001890528.vsix > jupyter.vsix && \
     bsdtar -xvf jupyter.vsix extension && \
     mv extension ms-toolsai.jupyter-2022.7.1001890528


RUN cd ${HOME}/.openvscode-server/extensions && \
    curl -L https://open-vsx.org/api/ms-toolsai/jupyter-renderers/1.0.9/file/ms-toolsai.jupyter-renderers-1.0.9.vsix  > jupyter-ren.vsix && \
     bsdtar -xvf jupyter-ren.vsix extension && \
     mv extension ms-toolsai.jupyter-renderers-1.0.9

# Configurar Jupyter widgets para VS Code
RUN mkdir -p /home/$USER/.openvscode-server/data/Machine \
    && echo '{\n\
    "jupyter.widgetScriptSources": [\n "jsdelivr.com", \n "unpkg.com"\n ],\n\
}' > /home/$USER/.openvscode-server/data/Machine/settings.json 


# Instalación de PyMFEM
RUN pip install --upgrade pip && \
    pip install --no-cache-dir ipykernel jupyter notebook ipywidgets mpi4py

RUN git clone https://github.com/mfem/PyMFEM.git 

RUN cd ${HOME}/PyMFEM && pip install ./ -C"with-parallel=Yes" 

RUN pip install matplotlib

# Instalación de websockets
RUN cd ${HOME} && curl -L https://github.com/aaugustin/websockets/archive/refs/tags/10.3.tar.gz > ./websockets-10.3.tar.gz && \
    tar xzf websockets-10.3.tar.gz && \
    cd websockets-10.3 && \
    pip install --user .

# GLVis JS
RUN git lfs install && git clone --depth=1 https://github.com/GLVis/glvis-js.git

RUN pip install glvis



# Puerto para VSCode Server (opcional)
EXPOSE 3000

USER root
ADD supervisord.conf /etc/supervisord.conf
RUN touch /var/log/glvis-js-webserver.log /var/log/glvis-browser-server.log /var/log/openvscode-server.log && \
    chown -R euler:euler /var/log/glvis-js-webserver.log /var/log/glvis-browser-server.log /var/log/openvscode-server.log

CMD ["/usr/bin/supervisord"]

