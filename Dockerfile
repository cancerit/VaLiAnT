FROM python:3.11-slim

# Install ps
RUN apt-get update && \
    apt-get install -yq --no-install-recommends procps

WORKDIR /tmp

# Install valiant package
COPY pyproject.toml setup.cfg .
COPY src/ src/
RUN pip install --no-cache-dir . && \
    rm -rf pyproject.toml setup.cfg src/

WORKDIR /home

CMD ["/bin/bash"]
