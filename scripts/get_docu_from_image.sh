#!/bin/bash -e
# script which copies files from the dockerimage $1(magnum.af.cpu-dev) into folders html and latex

if [ -n "$1" ]; then
    image=$1
    echo "using provided image $image"
else
    image=magnum.af.cpu-dev
    echo "using default image $image"
fi

id=$(docker create "$image")
docker cp $id:/home/magnum.af/docu/html/ - > html.tar
docker cp $id:/home/magnum.af/docu/latex/ - > latex.tar
docker rm -v $id

tar -xf html.tar
rm html.tar
echo "extracted $PWD/html/"

tar -xf latex.tar
rm latex.tar
echo "extracted $PWD/latex/"
