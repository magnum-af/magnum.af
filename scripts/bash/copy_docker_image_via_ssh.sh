#!/bin/bash -e
# script sending the docker image $1 to host $2 via ssh
# $1 image
# $2 user@host
# requirements: bzip2, pv, ssh
# apapted from https://stackoverflow.com/questions/23935141/how-to-copy-docker-images-from-one-host-to-another-without-using-a-repository
# example: ./$0 magnum.af user@host

if [ -n "$1" ]; then
    image=$1
    echo "using provided image $image"
else
    image=hello-world
    echo "using default image $image"
fi

if [ -n "$2" ]; then
    user_at_host=$2
    echo "using provided user_at_host $user_at_host"
else
    user_at_host=heistracher@131.130.24.70
    echo "using default user_at_host $user_at_host"
fi

docker save "$image" | bzip2 | pv | \
     ssh "$user_at_host" 'bunzip2 | docker load'
