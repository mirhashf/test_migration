#!/bin/bash

set -e

if ! pwd | grep "seqalto\/loomis2\/sh$"; then
    echo "This script should be run while you are in the seqalto/loomis2/sh directory."
    exit 1
fi

BINA_VER=$(./version-from-git)

./build-release

./release-package machine-binabox-backend default > bina_v${BINA_VER}.tar.gz
sha1sum bina_v${BINA_VER}.tar.gz > bina_v${BINA_VER}.tar.gz.SHA1SUM

./release-package portal-database default > bina-portal-database_v${BINA_VER}.tar.gz
sha1sum bina-portal-database_v${BINA_VER}.tar.gz > bina-portal-database_v${BINA_VER}.tar.gz.SHA1SUM

~/lake/opt/s3cmd/current/s3cmd -c ~/.s3cfg.DevEng.HenryChen --server-side-encryption put bina_v${BINA_VER}.tar.gz s3://bina.aws.dev-eng.deployment-archive/
~/lake/opt/s3cmd/current/s3cmd -c ~/.s3cfg.DevEng.HenryChen --server-side-encryption put bina_v${BINA_VER}.tar.gz.SHA1SUM s3://bina.aws.dev-eng.deployment-archive/

~/lake/opt/s3cmd/current/s3cmd -c ~/.s3cfg.DevEng.HenryChen --server-side-encryption put bina-portal-database_v${BINA_VER}.tar.gz s3://bina.aws.dev-eng.deployment-archive/
~/lake/opt/s3cmd/current/s3cmd -c ~/.s3cfg.DevEng.HenryChen --server-side-encryption put bina-portal-database_v${BINA_VER}.tar.gz.SHA1SUM s3://bina.aws.dev-eng.deployment-archive/

pushd ../../portal/templates

./generate
./package
tar -zcvf templates-${BINA_VER}.tar.gz templates-${BINA_VER}
sha1sum templates-${BINA_VER}.tar.gz > templates-${BINA_VER}.tar.gz.SHA1SUM

~/lake/opt/s3cmd/current/s3cmd -c ~/.s3cfg.DevEng.HenryChen --server-side-encryption put templates-${BINA_VER}.tar.gz s3://bina.aws.dev-eng.deployment-archive/templates/
~/lake/opt/s3cmd/current/s3cmd -c ~/.s3cfg.DevEng.HenryChen --server-side-encryption put templates-${BINA_VER}.tar.gz.SHA1SUM s3://bina.aws.dev-eng.deployment-archive/templates/

popd

