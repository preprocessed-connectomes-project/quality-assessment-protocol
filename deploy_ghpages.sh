#! /bin/bash
cd ./docs
jekyll build
mv ./_site ~/site
#cd ~ && git clone -b gh-pages ${CIRCLE_PROJECT_REPONAME} currentsite
cd ~ && git clone -b gh-pages_test ${CIRCLE_PROJECT_REPONAME} currentsite
mv ~/site/* ~/currentsite && cd ~/currentsite
git add -A
git commit -m "Automatic documentation update for build "${CIRCLE_SHA1}
git push
