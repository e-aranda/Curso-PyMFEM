# Instrucciones para construir una imagen Docker con PyMFEM
## Linux / MacOS
Es necesario tener instalado Docker previamente. Desde una terminal ejecutar:

```
$ sh docker-pymfem.sh download
$ sh docker-pymfem.sh run
```

O si se quiere construir la imagen:
```
sh docker-pymfem.sh build
```

Para ordenadores con procesador Apple Silicon usar
```
sh docker-pymfe-apple.sh download
sh docker-pymfe-apple.sh run
```
y para construir la imagen:

```
sh docker-pymfem-apple.sh build
```

## Windows
Es necesario tener instalado Docker previamente. Desde una terminal ejecutar:

```
$ ./docker-pymfem.ps1 download
$ ./docker-pymfem.ps1 run
```

O si se quiere construir la imagen:
```
$ ./docker-pymfem.ps1 build
```
