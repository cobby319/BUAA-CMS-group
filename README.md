Welcome to CMS group of BUAA
=======================================================================

There are some analysis framework and common codes

Instructions:
-----------------------
To clone this repository, you can use
```
git clone https://github.com/cobby319/BUAA-CMS-group.git
```
And then, to create a new branch
```
git checkout -b <branchname>
```
Add this repository by
```
git remote add <reponame you like to use> git@github.com:cobby319/BUAA-CMS-group.git
```
And some commands to commit your work
```
git add <file>
git commit -m "your comments"
```
Finally, push it to your branch
```
git push <reponame> <branchname>
```
B Physics test
=======================================================================

The recipe for analysis

```
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
cmsenv

git clone https://github.com/cobby319/BUAA-CMS-group.git -b b-phy-test
scram b -j 12

cd BUAA-CMS-group/JPsiUpsilonPhi/test/jpsipipi
vi test2016.py # you need change Line.24 input = cms.untracked.int32(-1), like 5000
cmsRun test2016.py
```
