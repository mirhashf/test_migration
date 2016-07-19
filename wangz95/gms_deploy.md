# GMS Deployment

## Prerequisites

Create a local clone of your fork of the `BinaTechnologies/seqalto` repo, then add `BinaTechnologies/seqalto` as remote called `upstream`
```bash
cd seqalto

git add remote upstream git@github.com:BinaTechnologies/seqalto.git
```

Fetch latest changes from `upstream` and merge `upstream/develop`
```bash
cd seqalto && git checkout develop && git fetch upstream && git merge upstream/develop && git push origin
```

Create a test tag
```bash
git tag -am"test-ctDNA-20160719" test-ctDNA-20160719

git push origin test-ctDNA-20160719

git push upstream test-ctDNA-20160719
```

SSH to `build-00` from `t-rex`
```bash
ssh build-00 -A
```

Clone the `seqalto` repository to your home directory on `build-00`, checkout a specific test tag:
```bash
git clone https://github.com/BinaTechnologies/seqalto.git

cd seqalto

git checkout test-ctDNA-20160719
```

Download 3rd-party packages using the Python script:
```bash
./third-party.py
```

## Installing

### In a shared environment

Set up a Python virtual environment
```bash
virtualenv venv

. ~/venv/bin/activate
```
Your shell prompt should now be prefixed with `(venv)`.
To return to the system Python you can run `deactivate`.

```bash
cd seqalto/tools/binaops
```

Run `setup.py`

```bash
python setup.py develop
```
After this you should have the utility called `deploy-gms`.

## GMS Build and Deploy

### Building and Packaging

Build and package using`deploy-gms`.

```bash
deploy-gms build
```

```bash
deploy-gms package
```
Make a note of the temporary directory returned from the packaging process.

### Adding the new packages to the APT repo

SSH to APT repo server `util-01` from `t-rex`
```bash
ssh util-01 -A
```

Copy `/home/troyerm/qe-repo-update.sh` to your own home directory on `util-01`:
```bash
cp -p /home/troyerm/qe-repo-update.sh ~/
```

Run `qe-repo-update.sh` with the temporary directory created by `deploy-gms package`:
```bash
./qe-repo-update.sh /tmp/tmpXXXXXX
```

If packages are not published correctly, run following command again:
```bash
rm -f debs/*

sudo -u binaops aptly publish update bina-qe bina
```

### Deploying

#### Pegging the GMS package version

Edit `/etc/salt/grains` on the target server
```bash
sudo vi /etc/salt/grains
```

Add an entry for `gms_iteration`
```
gms_iteration: 2977+g271a830
```
The value for `gms_iteration` is from the full version of packages built from `deploy-gms package`.

#### kent

SSH to `kent` from `t-rex`
```bash
ssh kent -A
```

Update the installed packages
```bash
sudo apt-get update && sudo salt-call state.highstate
```

#### tehran-20

SSH to `tehran-20` from `t-rex`
```bash
ssh tehran-20 -A
```

Update the installed packages
```bash
sudo apt-get update && sudo salt-call state.highstate
```

#### tehran-22

SSH to `tehran-22` from `t-rex`
```bash
ssh tehran-22 -A
```

Update the installed packages
```bash
sudo apt-get update && sudo salt-call state.highstate
```

Restart processes
```bash
sudo systemctl restart portal-frontend
sudo systemctl restart portal-backend
sudo systemctl restart executor
```

## References

[BinaOps Utilites](https://github.com/BinaTechnologies/seqalto/blob/develop/doc/gms/tools/BinaOps/BinaOps_utils.md)

[GMS Build and Deploy](https://github.com/BinaTechnologies/seqalto/blob/develop/doc/gms/tools/BinaOps/Build_and_Deploy.md)

