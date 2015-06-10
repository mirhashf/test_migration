#!/bin/bash

# TODO: pre-fill portal info if don't want prompt
PORTAL_URL=
PORTAL_USER=

if [ -z ${PORTAL_URL} ]; then
    read -p "Enter portal URL: " PORTAL_URL
fi

if [ -z ${PORTAL_USER} ]; then
    read -p "Enter portal username: " PORTAL_USER
fi

echo "Using portal ${PORTAL_URL} with user ${PORTAL_USER}"

CLI_JAR=$(mktemp)
UPDATE_SCRIPT=$(mktemp)

curl -fk ${PORTAL_URL}/cli/binaclient.jar > ${CLI_JAR} 2> /dev/null || FAIL_CLI_DL=1
echo 'jobobject.workflow.debug.keepScratch = "false"' > ${UPDATE_SCRIPT}

if [ ! -z "${FAIL_CLI_DL}" ]; then
    rm ${CLI_JAR}
    rm ${UPDATE_SCRIPT}

    echo "Error downloading CLI jar. Please check your portal URL."
    exit 1
fi

java -jar ${CLI_JAR} ${PORTAL_URL} ${PORTAL_USER} upgrade --template-js ${UPDATE_SCRIPT}

rm ${CLI_JAR}
rm ${UPDATE_SCRIPT}
