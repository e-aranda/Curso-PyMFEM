# Script para construir y ejecutar el contenedor MFEM usando Docker CLI en Windows

# Si hay problemas con la ejecución ejecutar desde una terminal PoweShell como
# administrador
# Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

# Configuración de salida (Colores)
$Green = "Green"
$Yellow = "Yellow"
$Red = "Red"

# Variables de configuración
$IMAGE_NAME = "earanda/pymfem:tutorial"
$CONTAINER_NAME = "pym"
$USER_NAME = "euler"
$USER_ID = 1000
$GROUP_ID = 1000

Write-Host "=== Script de construcción PyMFEM Docker (PowerShell) ===" -ForegroundColor $Green

# Función para mostrar uso
function Show-Usage {
    Write-Host "`nUso: .\build-and-run.ps1 [opcion]"
    Write-Host "`nOpciones:"
    Write-Host "  build      - Construir la imagen Docker"
    Write-Host "  run        - Ejecutar el contenedor"
    Write-Host "  exec       - Iniciar/Adjuntar al contenedor"
    Write-Host "  status     - Ver estado del contenedor"
    Write-Host "  download   - Descargar imagen oficial"
    exit
}

# Función para construir
function Build-Image {
    Write-Host "Construyendo imagen Docker..." -ForegroundColor $Yellow
    docker build `
        --build-arg USER=$USER_NAME `
        --build-arg UID=$USER_ID `
        --build-arg GID=$GROUP_ID `
        -t "${IMAGE_NAME}:tutorial" .
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host " Imagen construida exitosamente" -ForegroundColor $Green
        Write-Host "Imagen: ${IMAGE_NAME}" -ForegroundColor $Yellow
    }
}

# Función para ejecutar
function Run-Container {
    Write-Host "Iniciando contenedor..." -ForegroundColor $Yellow
    
    # Crear directorio workspace si no existe
    if (!(Test-Path "workspace")) {
        New-Item -ItemType Directory -Path "workspace" | Out-Null
    }

    # ${PWD} en PowerShell devuelve la ruta absoluta actual
    docker run --cap-add=SYS_PTRACE -p 3000:3000 -p 8000:8000 -p 8080:8080 `
        --name $CONTAINER_NAME `
        -v "${PWD}/workspace:/home/${USER_NAME}/workspace" `
        "${IMAGE_NAME}"
}

# Función para ejecutar en modo interactivo
function Exec-Container {
    Write-Host "Iniciando contenedor existente..." -ForegroundColor $Yellow
    docker start $CONTAINER_NAME
    docker attach $CONTAINER_NAME
}

# Función para ver estado
function Show-Status {
    Write-Host "Estado del contenedor:" -ForegroundColor $Yellow
    $container = docker ps -a --filter "name=$CONTAINER_NAME" --format "{{.Status}}"
    if ($container) {
        docker ps -a --filter "name=$CONTAINER_NAME"
    } else {
        Write-Host "Contenedor no encontrado" -ForegroundColor $Red
    }

    Write-Host "  Imágenes disponibles:" -ForegroundColor $Yellow
    docker images $IMAGE_NAME
}

# Función para descargar
function Download-Image {
    Write-Host "Descargando imagen de PyMFEM..." -ForegroundColor $Yellow
    docker pull earanda/pymfem:tutorial
}

# Procesar argumentos
if ($args.Count -eq 0) { Show-Usage }

switch ($args[0]) {
    "build"    { Build-Image }
    "run"      { Run-Container }
    "exec"     { Exec-Container }
    "status"   { Show-Status }
    "download" { Download-Image }
    Default    { Show-Usage }
}
