#!/bin/sh -ex
# adapted from https://pythonspeed.com/articles/faster-multi-stage-builds/

#docker build -t magnumaf_cpu_test -f Dockerfile.cpu . && ./bash/magnum.af.docker -t magnumaf_cpu_test python/examples/sp4.py -o output_sp4_docker_test

builder_image=magnumaf_cpu_test_builder
# docker pull "$builder_image" || true
docker build --target builder \
    --cache-from "$builder_image" \
    --tag "$builder_image" \
    --file Dockerfile.cpu .
# docker push "$builder_image"

runtime_image=magnumaf_cpu_test_runtime
# docker pull "$runtime_image" || true
docker build --target runtime \
    --cache-from "$runtime_image" \
    --tag "$runtime_image" \
    --file Dockerfile.cpu .
# docker push "$runtime_image"

./bash/magnum.af.docker --tag "$runtime_image" python/examples/sp4.py -o output_sp4_docker_test
