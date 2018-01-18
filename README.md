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
cmsrel CMSSW_9_2_3_patch2
cd CMSSW_9_2_3_patch2/src
cmsenv

git clone https://github.com/cobby319/BUAA-CMS-group.git
scram b -j 12

cd BUAA-CMS-group/JPsiUpsilonPhi/test
cmsRun miniAODmuonsRootupler.py
```
