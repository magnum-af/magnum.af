id=$(docker create magnum.af.cpu-dev)
docker cp $id:/home/magnum.af/html/ - > magnum.af_documentation_html.tar
docker cp $id:/home/magnum.af/latex/ - > magnum.af_documentation_latex.tar
docker rm -v $id
