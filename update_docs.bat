
REM Make sure we're on master branch
git checkout master

REM build docs
cd docs
call make clean
call make html
cd ..

REM commit and push
git add .
git commit
git push origin master

REM switch branches and pull the data we want
git checkout gh-pages
git rm *.rst
git rm *.html
git rm *.inv
git rm *.js
git rm *.buildinfo
git rm .\_sources\*
git rm .\_static\*
git rm .\_modules\*
git rm docs -rf

REM Move it to the root and remove the docs directory
git checkout master docs/build/html
call xcopy .\docs\build\html .\ /E
git rm docs -rf

REM Add and push
git add .
git commit
git push origin gh-pages

REM Switch back to master
git checkout master