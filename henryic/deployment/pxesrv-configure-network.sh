#!/bin/bash

IFACE=eth1

NETMASK=255.255.255.0
CIDR=24
GATEWAY_IP=10.0.0.1
DNS_IPS=10.0.0.1
DOMAIN=binatechnologies.com

function usage {
    echo "usage: ${0} <up|down>"
    exit 1
}

function bring_up {
    echo "Bringing up"
    ip link set dev ${IFACE} up
    ip address add ${GATEWAY_IP}/${CIDR} dev ${IFACE}

    LABEL_COUNT=0
    for IP in $(echo ${DNS_IPS} | tr "," "\n"); do
        ip address add ${IP}/${CIDR} dev ${IFACE} label ${IFACE}:${LABEL_COUNT}
        ((++LABEL_COUNT))
    done
}

function tear_down {
    echo "Tearing down"

    for IP in $(ip address show eth1 | grep "^ \+inet " | grep "${IFACE}:" | awk '{ print $2 }'); do
        ip address del ${IP} dev ${IFACE}
    done

    ip address del ${GATEWAY_IP}/${CIDR} dev ${IFACE}

    ip link set dev ${IFACE} down
}

[ ! -z ${1} ] || usage

if [ "${1}" == "up" ]; then
    bring_up
elif [ "${1}" == "down" ]; then
    tear_down
else
    usage
fi

