#! /bin/bash
#cd ./docs
#jekyll build
#mv ./_site ~/site
#cd ~ && git clone -b gh-pages ${CIRCLE_PROJECT_REPONAME} currentsite
cd ~ && git clone -b gh-pages_test https://github.com/preprocessed-connectomes-project/quality-assessment-protocol.git currentsite
#mv ~/site/* ~/currentsite && cd ~/currentsite
cd ~/currensite
touch ${CIRCLE_SHA1}
git add -A
git commit -m "Automatic documentation update for build "${CIRCLE_SHA1}
git push
