# Imagen de loxofed. Construir desde la raiz del repositorio:
#   docker build -f docker/Dockerfile -t loxofed:0.4.0 .
# rocker/r-ver fija el snapshot de CRAN a la fecha de la imagen: mismas versiones de paquetes
# en los 7 sitios. NO se instala nada en tiempo de ejecucion (el contenedor corre sin red).
FROM rocker/r-ver:4.4.2
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

RUN install2.r --error --skipinstalled data.table bnlearn yaml jsonlite digest testthat \
    && rm -rf /tmp/downloaded_packages

COPY . /src/loxofed
# --install-tests deja los tests dentro de la imagen para poder verificarla:
#   docker run --rm --entrypoint Rscript loxofed:0.4.0 -e 'testthat::test_package("loxofed")'
RUN R CMD INSTALL --install-tests /src/loxofed && rm -rf /src

# Usuario sin privilegios; las tres rutas de trabajo se montan desde el anfitrion
RUN useradd -m -u 1000 loxo && mkdir -p /data /config /output && chown loxo /output
USER loxo
ENV HOME=/tmp
WORKDIR /output
VOLUME ["/data", "/config", "/output"]

ENTRYPOINT ["Rscript", "-e", "loxofed::lox_cli()"]
CMD ["help"]
