#!/bin/bash -e
id=$(docker create magnum.af.cpu-dev)
docker cp $id:/home/magnum.af/docu/html/ - > html.tar
docker cp $id:/home/magnum.af/docu/latex/ - > latex.tar
docker rm -v $id

tar -xf html.tar
rm html.tar
echo "extracted $PWD/html/"

tar -xf latex.tar
rm latex.tar
echo "extracted $PWD/latex/"
