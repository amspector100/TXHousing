
# Make sure we're on master branch
git checkout master

# build docs
cd docs
make clean
make html
cd ..

# commit and push
git add .
git commit
git push origin master

# switch branches and pull the data we want
git checkout gh-pages
git rm .rst
git rm .html
git rm docs -rf

# Move it to the root and remove the docs directory
git checkout master docs/build/html
xcopy .\docs\build\html .\
git rm docs -rf

# Add and push
git add .
git commit
git push origin gh-pages

# Switch back to master
git checkout master