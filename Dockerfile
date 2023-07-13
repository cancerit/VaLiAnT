FROM python:3.9-slim AS builder

ENV CYTHON_VERSION 0.29.30

# Install build dependencies
RUN apt-get update && \
    apt-get install -yq --no-install-recommends build-essential && \
    pip install --no-cache-dir cython==${CYTHON_VERSION}

WORKDIR /tmp

# Build Python package wheels
WORKDIR /wheels
COPY src/requirements.txt .
RUN pip wheel --no-cache-dir -r requirements.txt && \
    rm requirements.txt

FROM python:3.9-slim

# Install ps
RUN apt-get update && \
    apt-get install -yq --no-install-recommends procps

WORKDIR /tmp

# Install wheels
COPY --from=builder /wheels /wheels
RUN pip install --no-warn-script-location /wheels/*.whl && \
    rm -rf /wheels && \
    rm -rf /root/.cache/pip/*

# Install valiant package
COPY src/ src/
RUN pip install --no-cache-dir src/ && \
    rm -rf src/

WORKDIR /home

CMD ["/bin/bash"]
