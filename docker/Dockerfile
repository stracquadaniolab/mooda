FROM  python:3.7-slim

LABEL org.opencontainers.image.title="mooda"
LABEL org.opencontainers.image.description="Multi-Objective Optimization for DNA design and assembly"
LABEL org.opencontainers.image.url="https://github.com/stracquadaniolab/mooda"
LABEL org.opencontainers.image.documentation="https://github.com/stracquadaniolab/mooda"
LABEL org.opencontainers.image.source="https://github.com/stracquadaniolab/mooda"
LABEL org.opencontainers.image.revision="v0.9.2"
LABEL org.opencontainers.image.vendor="stracquadaniolab"
LABEL org.opencontainers.image.authors="Giovanni Stracquadanio"

COPY . /app

# Install basic tools
RUN apt-get update \
    && apt-get install --yes $(cat app/docker/pkg.list) \
    && apt-get autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN cd /app && python setup.py install


ENTRYPOINT ["tini", "--", "/usr/local/bin/mooda"]
