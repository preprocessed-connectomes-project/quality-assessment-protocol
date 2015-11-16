#! /bin/bash
#cd ./docs
#jekyll build
mv ./_site ~/site
cd ~ && git clone -b gh-pages_test https://github.com/preprocessed-connectomes-project/quality-assessment-protocol.git currentsite
mv ~/site/* ~/currentsite && cd ~/currentsite
git add -A
git commit -m "Automatic documentation update for build "${CIRCLE_SHA1}
git push
