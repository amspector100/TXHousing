# build docs
cd docs
make clean
make html
cd ..

# commit and push
git add .
git commit -m "building and pushing docs"
git push origin master

# switch branches and pull the data we want
git checkout gh-pages
git rm -rf . # (Clear gh-pages branch)
touch .nojekyll
git checkout master docs/build/html

# Move it to the project root and also get rid of the docs file
xcopy '.\docs\build\html' '.\'
git rm -rf docs\.

git add .
git commit -m "publishing updated docs"
git push origin gh-pages
