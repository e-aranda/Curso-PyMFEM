#!/bin/bash

# Script para construir y ejecutar el contenedor MFEM usando Docker CLI

set -e

# Colores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Variables de configuración
IMAGE_NAME="earanda/pymfem:tutorial
"
CONTAINER_NAME="pym"
USER_NAME="euler"
USER_ID=1000
GROUP_ID=1000

echo -e "${GREEN}=== Script de construcción PyMFEM Docker (CLI) ===${NC}\n"

# Función para mostrar uso
usage() {
    echo "Uso: $0 [opcion]"
    echo ""
    echo "Opciones:"
    echo "  build      - Construir la imagen Docker"
    echo "  run        - Ejecutar el contenedor (detached)"
    echo "  shell      - Abrir shell en el contenedor"
    echo "  status     - Ver estado del contenedor"
    echo "  download - Descargar imagen"
    exit 1
}

# Función para construir
build_image() {
    echo -e "${YELLOW}Construyendo imagen Docker...${NC}"
    docker build \
        --build-arg USER=${USER_NAME} \
        --build-arg UID=${USER_ID} \
        --build-arg GID=${GROUP_ID} \
        -t ${IMAGE_NAME} \
        .
    echo -e "${GREEN}✓ Imagen construida exitosamente${NC}"
    echo -e "${YELLOW}Imagen: ${IMAGE_NAME}${NC}"
}

# Función para ejecutar en modo detached
run_container() {
    echo -e "${YELLOW}Iniciando contenedor en modo detached...${NC}"
    
    # Crear directorios si no existen
    mkdir -p workspace
    docker run --cap-add=SYS_PTRACE -p 3000:3000 -p 8000:8000 -p 8080:8080 \
	   --name ${CONTAINER_NAME} \
	    -v "$(pwd)/workspace:/home/${USER_NAME}/workspace" ${IMAGE_NAME}

    
}

# Función para ejecutar en modo interactivo
exec_container() {
    echo -e "${YELLOW}Iniciando contenedor...${NC}"
    
    # Crear directorios si no existen
    docker start ${CONTAINER_NAME}
    docker attach ${CONTAINER_NAME}
}

# Función para ver estado
show_status() {
    echo -e "${YELLOW}Estado del contenedor:${NC}"
    docker ps -a | grep ${CONTAINER_NAME} || echo -e "${RED}Contenedor no encontrado${NC}"
    echo ""
    echo -e "${YELLOW}Imágenes disponibles:${NC}"
    docker images | grep ${IMAGE_NAME} || echo -e "${RED}Imagen no encontrada${NC}"
}

# Función para usar imagen oficial
download_image() {
    echo -e "${YELLOW}Descargando imagen de PyMFEM...${NC}"
    docker pull earanda/pymfem:tutorial
}



# Procesar argumentos
case "$1" in
    build)
        build_image
        ;;
    run)
        run_container
        ;;
    exec)
        exec_container
        ;;
    status)
        show_status
        ;;
    download)
	download_image
	;;
    *)
        usage
        ;;
esac
