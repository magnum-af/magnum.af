#!/bin/bash -e
# retags docker image $1 with $2 when [ "$3" == "master" ]

if [ "$3" == "master" ]; then
    echo "on branch master: performing pull, retag, push"
    docker pull $1
    docker tag $1 $2
    docker push $2
else
    echo "not on master, no retag"
fi
