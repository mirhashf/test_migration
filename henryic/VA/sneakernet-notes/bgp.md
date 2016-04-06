Setting up BGP on VA Transfer Server
====================================

- Ensure entire link chain is set to auto-negotiate
- Install quagga + dependencies
- Install ethtool
- Install vlan

- Set eth1 to VLAN 84
- Set ports on switch to VLAN 84


Configure quagga
----------------

/etc/quagga/bgpd.conf:

>    ! -*- bgp -*-
>    !
>    ! BGPd sample configuratin file
>    !
>    ! $Id: bgpd.conf.sample,v 1.1 2002/12/13 20:15:29 paul Exp $
>    !
>    hostname bgpd
>    password zebra
>    !enable password please-set-at-here
>    !
>    !bgp mulitple-instance
>    !
>    router bgp 65448
>      neighbor 67.223.179.6 remote-as 25658
>      neighbor 67.223.179.6 password <BGP MD5 password>
>    ! bgp router-id 10.0.0.1
>    ! network 10.0.0.0/8
>    ! neighbor 10.0.0.2 remote-as 7675
>    ! neighbor 10.0.0.2 route-map set-nexthop out
>    ! neighbor 10.0.0.2 ebgp-multihop
>    ! neighbor 10.0.0.2 next-hop-self
>    !
>    ! access-list all permit any
>    !
>    !route-map set-nexthop permit 10
>    ! match ip address all
>    ! set ip next-hop 10.0.0.1
>    !
>    !log file /var/log/quagga/bgpd.log
>    !
>    log stdout

/etc/quagga/daemons:
>    # This file tells the quagga package which daemons to start.
>    #
>    # Entries are in the format: <daemon>=(yes|no|priority)
>    #   0, "no"  = disabled
>    #   1, "yes" = highest priority
>    #   2 .. 10  = lower priorities
>    # Read /usr/share/doc/quagga/README.Debian for details.
>    #
>    # Sample configurations for these daemons can be found in
>    # /usr/share/doc/quagga/examples/.
>    #
>    # ATTENTION:
>    #
>    # When activation a daemon at the first time, a config file, even if it is
>    # empty, has to be present *and* be owned by the user and group "quagga", else
>    # the daemon will not be started by /etc/init.d/quagga. The permissions should
>    # be u=rw,g=r,o=.
>    # When using "vtysh" such a config file is also needed. It should be owned by
>    # group "quaggavty" and set to ug=rw,o= though. Check /etc/pam.d/quagga, too.
>    #
>    # The watchquagga daemon is always started. Per default in monitoring-only but
>    # that can be changed via /etc/quagga/debian.conf.
>    #
>    zebra=yes
>    bgpd=yes
>    ospfd=no
>    ospf6d=no
>    ripd=no
>    ripngd=no
>    isisd=no
>    babeld=no

/etc/quagga/vtysh.conf:
>    !
>    ! Sample
>    !
>    ! service integrated-vtysh-config
>    hostname sneaker
>    username root nopassword
>    !

